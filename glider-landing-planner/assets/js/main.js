import { clampHandle, validHandle } from './model/limits.js?v=0.1.20';
import { isCompactLayout, COMPACT_LAYOUT_QUERY } from './ui/responsive.js?v=0.1.20';
import { createDemo, makeStraight, makeTurn, trafficPreset, ENGINE_VERSION } from './model/defaults.js?v=0.1.20';
import { buildPlan, rememberResolvedTurnAngles } from './planning/plan.js?v=0.1.20';
import { chooseConnector } from './planning/connector.js?v=0.1.20';
import { validateProject, loadProject, parseProject, downloadProject } from './io/project.js?v=0.1.20';
import { saveSnapshot, latestSnapshot } from './io/autosave.js?v=0.1.20';
import { PlannerMap } from './ui/map.js?v=0.1.20';
import { chartSvg, fmtAlt, fmtElapsed } from './ui/chart.js?v=0.1.20';
import { openPrintReport } from './ui/print.js?v=0.1.20';
import { checkpointTable } from './ui/checkpoints.js?v=0.1.20';
import { heightInputText, FEET_PER_METER } from './ui/height-units.js?v=0.1.20';
import { escapeHtml, htmlLegId, readLegId } from './ui/escape-html.js?v=0.1.20';
import { effectiveLegSpeedMps, effectiveStraightDirectionDeg, effectiveTurnValues } from './ui/leg-display.js?v=0.1.20';
const n = (value, d = 2) => Number.isFinite(value) ? Number(value.toFixed(d)) : 0;
const options = (items, selected) => items.map(([v, t]) => `<option value="${v}" ${v === selected ? 'selected' : ''}>${t}</option>`).join('');
const engineNote = '教育・机上検討用。実際の飛行可否、地形・障害物、フレア・スポイラー性能は評価しません。';
export class App {
    root;
    project = createDemo();
    plan = buildPlan(this.project);
    history = [];
    future = [];
    map;
    selected = null;
    pick = null;
    saveTimer;
    status = '未保存';
    notice = '';
    editRevision = 0;
    mobilePanel = null;
    constructor(root) {
        this.root = root;
        root.innerHTML = `<div class="glp-shell"><div class="glp-toolbar"></div><div class="glp-app-body"><aside class="glp-sidebar"><div class="glp-sidebar-scroll"></div></aside><section class="glp-main"><div class="glp-mapbar"><div class="glp-map-actions"><button data-action="pick-departure" class="glp-btn glp-btn-quiet">離脱地点を指定</button><button data-action="pick-landing" class="glp-btn glp-btn-quiet">着陸地点を指定</button><button data-action="fit" class="glp-btn glp-btn-quiet">全体表示</button></div><span class="glp-map-hint">地図上の経路・地点を選択して編集できます。紫のBezier制御ハンドルのみドラッグ可能です</span></div><div class="glp-map"></div><div class="glp-overview"><div class="glp-overview-top"><strong>高度プロファイル</strong><span class="glp-chart-legend"><i class="glp-line-blue"></i>離脱側予想高度 <i class="glp-line-green"></i>着陸側逆算必要高度</span></div><div class="glp-chart"></div><div class="glp-metrics"></div><div class="glp-checkpoints"></div><div class="glp-diagnostics"></div></div></section></div><div class="glp-footnote">${engineNote}　｜　地図：<a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>（地理院タイル）</div><input class="glp-import-file" type="file" accept=".json,.glp.json,application/json" hidden><section class="glp-mobile-editor" hidden aria-label="地図から経路を編集"><div class="glp-mobile-editor-head"></div><div class="glp-mobile-editor-content"></div></section></div>`;
        this.map = new PlannerMap(root.querySelector('.glp-map'), {
            onPick: (type, lat, lon) => {
                this.pick = null;
                this.map.setPick(null);
                this.selected = null;
                this.edit(p => {
                    if (type === 'departure') {
                        p.departure.latDeg = lat;
                        p.departure.lonDeg = lon;
                    }
                    else {
                        p.landing.point.latDeg = lat;
                        p.landing.point.lonDeg = lon;
                    }
                });
                this.openMobilePanel({ kind: 'point', point: type });
                this.scrollSidebarToPoint(type);
            },
            onHandle: (end, lengthM) => {
                const changed = this.edit(p => {
                    if (p.connector) {
                        p.connector.adjustmentMode = 'manual';
                        if (end === 'start')
                            p.connector.handleStartM = clampHandle(lengthM);
                        else
                            p.connector.handleEndM = clampHandle(lengthM);
                    }
                });
                if (!changed)
                    this.map.redraw(this.project, this.plan, this.selected);
            },
            onSelectLeg: (id) => this.selectLeg(id, true),
            onSelectPoint: (point) => { this.selected = null; this.renderSidebar(); this.map.redraw(this.project, this.plan, null); this.updateCheckpointSelection(); this.openMobilePanel({ kind: 'point', point }); this.scrollSidebarToPoint(point); },
            onAddAtEnd: (side, type) => {
                if (type)
                    this.handleAction(type === 'straight' ? 'add-straight' : 'add-turn', { side });
                else
                    this.openMobilePanel({ kind: 'add', side });
            }
        });
        root.addEventListener('click', e => {
            const btn = e.target.closest('button[data-action]');
            if (btn) {
                e.preventDefault();
                this.handleAction(btn.dataset.action ?? '', { ...btn.dataset, id: readLegId(btn.dataset.id) });
            }
        });
        root.addEventListener('change', e => {
            const el = e.target;
            if (el.matches('[data-path]'))
                this.handleChange(el);
            else if (el.classList.contains('glp-import-file'))
                void this.doImport(el);
        });
        this.render();
        this.map.showNewProject(this.project);
        // Redraw touch hit areas and endpoint + buttons when rotating or resizing a device.
        window.matchMedia?.(COMPACT_LAYOUT_QUERY).addEventListener?.('change', () => { this.map.redraw(this.project, this.plan, this.selected); this.renderMobilePanel(); });
        const revisionAtStartup = this.editRevision;
        void latestSnapshot().then(p => {
            if (!p || this.editRevision !== revisionAtStartup)
                return; // Don't overwrite edits made while IndexedDB was loading.
            let valid, plan;
            try {
                valid = loadProject(p);
                plan = this.prepareProject(valid);
            }
            catch (e) {
                this.notice = '自動保存データを復元できません（新規計画を表示）: ' + (e instanceof Error ? e.message : String(e));
                this.renderResults();
                return;
            }
            this.plan = plan;
            this.project = valid;
            this.history = [];
            this.future = [];
            this.selected = null;
            this.status = '自動保存から復元しました';
            this.mobilePanel = null;
            this.render();
            this.map.fit(this.project, this.plan);
        }).catch(() => { });
    }
    getField(path) {
        let a = this.project;
        for (const k of path.split('.'))
            a = a[k];
        return a;
    }
    control(path, title, unit = '', step = 'any', scale = 1, shownWhenInherited) {
        const raw = this.getField(path);
        const isHeight = path === 'departure.altitudeMSLM' || path === 'landing.elevationMSLM';
        const heightUnit = this.project.display.altitudeUnit;
        const actualScale = isHeight && heightUnit === 'ft' ? FEET_PER_METER : scale;
        const fieldUnit = isHeight ? `${heightUnit}, MSL` : unit;
        const displayed = (raw === null || raw === undefined) && shownWhenInherited !== undefined ? shownWhenInherited : raw;
        const value = isHeight && typeof displayed === 'number' ? heightInputText(displayed, heightUnit) : typeof displayed === 'number' ? n(displayed * scale, scale === 3.6 ? 2 : 5) : displayed ?? '';
        return `<label class="glp-field"><span>${title}${fieldUnit ? ` <em>${fieldUnit}</em>` : ''}</span><input type="number" step="${step}" value="${escapeHtml(value)}" data-path="${path}" data-scale="${actualScale}"></label>`;
    }
    checkbox(path, title) { return `<label class="glp-check"><input type="checkbox" data-path="${path}" ${this.getField(path) ? 'checked' : ''}><span>${title}</span></label>`; }
    opacityControl(path, title) { const value = Math.round(Number(this.getField(path)) * 100); return `<label class="glp-opacity-field"><span>${title} <output>${value}%</output></span><input type="range" min="0" max="100" step="5" value="${value}" data-path="${path}" data-scale="100"></label>`; }
    select(path, title, items) { return `<label class="glp-field"><span>${title}</span><select data-path="${path}">${options(items, String(this.getField(path)))}</select></label>`; }
    legEditor(side, leg, i) {
        const modeOpts = [['heading', '機首方位指定'], ['ground-track', '進行方向指定']];
        const base = `${side === 'departure' ? 'departureLegs' : 'arrivalLegs'}.${i}`;
        // Nullable model values retain their inherited meaning even when the form shows a number.
        const speedInForce = effectiveLegSpeedMps(this.project, leg);
        const directionInForce = effectiveStraightDirectionDeg(this.project, this.plan, side, leg);
        const turnValues = leg.type === 'turn' ? effectiveTurnValues(this.plan, side, leg) : null;
        const directionName = leg.guidance === 'heading' ? '機首方位' : '進行方向';
        const turnAngleLabel = `${directionName}変化量`;
        const turnTargetLabel = `旋回${side === 'departure' ? '後' : '前'}の${directionName}`;
        const turnAngleField = leg.type === 'turn' && leg.angleInputMode === 'target'
            ? `<label class="glp-field"><span>${turnAngleLabel} <em>°</em></span><input type="number" readonly value="${turnValues.completed ? n(turnValues.angleDeg, 5) : ''}" placeholder="評価不可"></label>`
            : this.control(base + '.angleDeg', turnAngleLabel, '°');
        const targetField = leg.type === 'turn' && leg.angleInputMode === 'target'
            ? this.control(base + '.targetDirectionDeg', turnTargetLabel, '°')
            : `<label class="glp-field"><span>${turnTargetLabel} <em>°</em></span><input type="number" readonly value="${turnValues?.boundaryDeg === null || turnValues === null ? '' : n(turnValues.boundaryDeg, 5)}" placeholder="評価不可"></label>`;
        return `<div class="glp-leg-edit"><div class="glp-two">${this.select(base + '.guidance', '飛行方式', modeOpts)}${this.control(base + '.speedOverrideMps', leg.type === 'turn' ? '旋回速度' : '滑空速度', 'km/h', '0.1', 3.6, speedInForce)}</div>
    ${leg.type === 'straight' ? `<div class="glp-two">${this.control(base + '.lengthM', '直線距離', 'm')}${this.control(base + '.directionDeg', directionName, '°', 'any', 1, directionInForce)}</div>` : `<div class="glp-two">${this.select(base + '.direction', '旋回方向', [['right', '右'], ['left', '左']])}${this.control(base + '.bankDeg', 'バンク角', '°')}</div>${this.select(base + '.angleInputMode', '入力方式', [['angle', `${turnAngleLabel}を指定`], ['target', `${directionName}を指定`]])}<div class="glp-two">${turnAngleField}${targetField}</div>`}
    <label class="glp-field"><span>チェックポイント名</span><input type="text" maxlength="70" value="${escapeHtml(leg.checkpointLabel ?? '')}" data-path="${base}.checkpointLabel"></label></div>`;
    }
    legList(side) {
        const legs = side === 'departure' ? this.project.departureLegs : this.project.arrivalLegs;
        const lines = legs.map((leg, i) => {
            const active = this.selected === leg.id;
            const effectiveAngle = leg.type === 'turn' ? effectiveTurnValues(this.plan, side, leg).angleDeg : 0;
            const title = leg.type === 'straight' ? `直線 ${Math.round(leg.lengthM)} m` : `${leg.direction === 'right' ? '右' : '左'}旋回 ${n(effectiveAngle, 2)}° / バンク${leg.bankDeg}°`;
            return `<div class="glp-leg ${active ? 'glp-leg-active' : ''}" data-leg-id="${htmlLegId(leg.id)}"><div class="glp-leg-line"><button class="glp-leg-title" data-action="select-leg" data-id="${htmlLegId(leg.id)}" title="編集"><span class="glp-leg-symbol">${leg.type === 'straight' ? '━' : '⤷'}</span><span>${i + 1}. ${escapeHtml(title)}</span></button><span class="glp-leg-buttons"><button title="上に移動" data-action="move-up" data-side="${side}" data-id="${htmlLegId(leg.id)}">↑</button><button title="下に移動" data-action="move-down" data-side="${side}" data-id="${htmlLegId(leg.id)}">↓</button><button title="削除" data-action="delete-leg" data-side="${side}" data-id="${htmlLegId(leg.id)}">×</button></span></div>
       ${active ? this.legEditor(side, leg, i) + `<div class="glp-desktop-leg-actions"><button class="glp-btn glp-btn-sub glp-desktop-remove" data-action="delete-leg" data-side="${side}" data-id="${htmlLegId(leg.id)}">このレグを削除</button></div>` : ''}</div>`;
        }).join('');
        const adds = `<div class="glp-adds"><button class="glp-btn glp-btn-sub" data-action="add-straight" data-side="${side}">＋ 直線</button><button class="glp-btn glp-btn-sub" data-action="add-turn" data-side="${side}">＋ 旋回</button></div>`;
        const list = `<div class="glp-leg-list">${lines || '<p class="glp-muted">まだレグがありません。</p>'}</div>`;
        return side === 'arrival' ? adds + list : list + adds;
    }
    selectLeg(id, fromMap = false) {
        const next = this.selected === id && !fromMap ? null : id;
        this.selected = next;
        this.renderSidebar();
        this.map.redraw(this.project, this.plan, this.selected);
        this.updateCheckpointSelection();
        if (next && isCompactLayout())
            this.openMobilePanel({ kind: 'leg', id: next });
        else if (next && fromMap)
            this.scrollSidebarTo(Array.from(this.root.querySelectorAll('.glp-leg[data-leg-id]')).find(el => readLegId(el.dataset.legId) === next) ?? null);
        else if (!next && this.mobilePanel?.kind === 'leg')
            this.openMobilePanel(null);
    }
    /** Scroll only the desktop editor, never the surrounding WordPress page. */
    scrollSidebarTo(target) {
        if (!target || isCompactLayout())
            return;
        const sidebar = this.root.querySelector('.glp-sidebar-scroll');
        if (!sidebar)
            return;
        const viewport = sidebar.getBoundingClientRect();
        const section = target.getBoundingClientRect();
        const padding = 12;
        if (section.top >= viewport.top + padding && section.bottom <= viewport.bottom - padding)
            return;
        const top = Math.max(0, sidebar.scrollTop + section.top - viewport.top - padding);
        const reduceMotion = window.matchMedia?.('(prefers-reduced-motion: reduce)')?.matches;
        if (typeof sidebar.scrollTo === 'function')
            sidebar.scrollTo({ top, behavior: reduceMotion ? 'auto' : 'smooth' });
        else
            sidebar.scrollTop = top;
    }
    scrollSidebarToPoint(point) {
        this.scrollSidebarTo(this.root.querySelector(`[data-point-editor="${point}"]`));
    }
    openMobilePanel(target) {
        if (!isCompactLayout())
            return;
        this.mobilePanel = target;
        const content = this.root.querySelector('.glp-mobile-editor-content');
        if (content)
            content.scrollTop = 0;
        this.renderMobilePanel();
    }
    renderMobilePanel() {
        const panel = this.root.querySelector('.glp-mobile-editor');
        if (!panel)
            return;
        if (!this.mobilePanel || !isCompactLayout()) {
            panel.hidden = true;
            return;
        }
        const target = this.mobilePanel;
        let title = '';
        let contents = '';
        if (target.kind === 'leg') {
            const departureIndex = this.project.departureLegs.findIndex(l => l.id === target.id);
            const side = departureIndex >= 0 ? 'departure' : 'arrival';
            const index = departureIndex >= 0 ? departureIndex : this.project.arrivalLegs.findIndex(l => l.id === target.id);
            if (index < 0) {
                this.mobilePanel = null;
                panel.hidden = true;
                return;
            }
            const leg = (side === 'departure' ? this.project.departureLegs : this.project.arrivalLegs)[index];
            title = `${side === 'departure' ? '離脱側' : '着陸側'} ${leg.type === 'turn' ? '旋回' : '直線'} ${index + 1}`;
            const list = side === 'departure' ? this.project.departureLegs : this.project.arrivalLegs;
            contents = this.legEditor(side, leg, index) + `<div class="glp-mobile-leg-actions"><button class="glp-btn glp-btn-quiet" data-action="move-up" data-side="${side}" data-id="${htmlLegId(leg.id)}" ${index === 0 ? 'disabled' : ''}>↑ 上に移動</button><button class="glp-btn glp-btn-quiet" data-action="move-down" data-side="${side}" data-id="${htmlLegId(leg.id)}" ${index === list.length - 1 ? 'disabled' : ''}>↓ 下に移動</button><button class="glp-btn glp-btn-sub glp-mobile-remove" data-action="delete-leg" data-side="${side}" data-id="${htmlLegId(leg.id)}">このレグを削除</button></div>`;
        }
        else if (target.kind === 'point') {
            if (target.point === 'departure') {
                title = '離脱地点・条件';
                contents = `<div class="glp-two">${this.control('departure.latDeg', '緯度', '°')}${this.control('departure.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('departure.altitudeMSLM', '離脱高度')}${this.control('departure.headingDeg', '開始機首方位', '°')}</div><button class="glp-btn glp-btn-sub" data-action="open-add" data-side="departure">＋ 離脱側にレグを追加</button>`;
            }
            else {
                title = '着陸地点・ファイナル';
                contents = `<div class="glp-two">${this.control('landing.point.latDeg', '緯度', '°')}${this.control('landing.point.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('landing.elevationMSLM', '着陸標高')}${this.control('landing.runwayTrackDeg', '真方位（進入）', '°')}</div><div class="glp-two">${this.control('landing.final.pathAngleDeg', '目標パス角', '°')}${this.control('landing.final.groundLengthM', '直線ファイナル長', 'm')}</div><p class="glp-tip">入力後、場周プリセットを追加できます（既存レグは確認のうえ置換）。</p><div class="glp-preset"><button class="glp-btn glp-btn-soft" data-action="preset-left">左場周プリセット</button><button class="glp-btn glp-btn-soft" data-action="preset-right">右場周プリセット</button></div><button class="glp-btn glp-btn-sub" data-action="open-add" data-side="arrival">＋ 着陸側にレグを追加</button>`;
            }
        }
        else {
            title = target.side === 'departure' ? '離脱側にレグを追加' : '着陸側にレグを追加';
            contents = `<p class="glp-tip">${target.side === 'departure' ? '離脱側の最後' : '着陸側の最上流'}に追加します。</p><div class="glp-adds"><button class="glp-btn glp-btn-sub" data-action="add-straight" data-side="${target.side}">＋ 直線</button><button class="glp-btn glp-btn-sub" data-action="add-turn" data-side="${target.side}">＋ 旋回</button></div>${target.side === 'departure' && !this.project.departureLegs.length ? '<button class="glp-btn glp-btn-quiet" data-action="edit-departure">離脱地点・条件を編集</button>' : ''}`;
        }
        panel.hidden = false;
        panel.querySelector('.glp-mobile-editor-head').innerHTML = `<strong>${escapeHtml(title)}</strong><button type="button" class="glp-btn glp-btn-quiet" data-action="close-mobile-panel" aria-label="編集画面を閉じる">閉じる ×</button>`;
        const body = panel.querySelector('.glp-mobile-editor-content');
        const lastScroll = body.scrollTop;
        body.innerHTML = contents;
        body.scrollTop = lastScroll;
    }
    renderSidebar() {
        const sidebarScroll = this.root.querySelector('.glp-sidebar-scroll');
        const preservedScroll = sidebarScroll.scrollTop;
        sidebarScroll.innerHTML = `
   <div class="glp-panel-heading"><div class="glp-eyebrow">FLIGHT PARAMETERS</div><h2>飛行条件・経路設計</h2></div>
   <div class="glp-panel"><h3>基本条件</h3><div class="glp-two">${this.control('defaultSpeedMps', '基準対気速度', 'km/h', '0.1', 3.6)}${this.control('maxBankDeg', '最大バンク', '°')}</div>
   <div class="glp-two">${this.control('performance.bestLD', '最大滑空比')}${this.control('performance.bestLDAtMps', '最良滑空速度', 'km/h', '0.1', 3.6)}</div>
   <div class="glp-two">${this.control('wind.fromDeg', '風向', '°')}${this.control('wind.speedMps', '風速', 'm/s')}</div>
   ${this.select('display.altitudeUnit', '高度入力・表示単位（MSL）', [['ft', 'ft'], ['m', 'm']])}</div>
   <div class="glp-panel" data-point-editor="departure"><h3><span class="glp-dot blue"></span>離脱地点・条件</h3><div class="glp-two">${this.control('departure.latDeg', '緯度', '°')}${this.control('departure.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('departure.altitudeMSLM', '離脱高度')}${this.control('departure.headingDeg', '開始機首方位', '°')}</div><p class="glp-tip">地図の「離脱地点を指定」で位置を選ぶか、緯度・経度を入力してください。</p></div>
   <div class="glp-panel" data-point-editor="landing"><h3><span class="glp-dot green"></span>着陸地点・ファイナル・場周設定</h3><div class="glp-two">${this.control('landing.point.latDeg', '緯度', '°')}${this.control('landing.point.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('landing.elevationMSLM', '着陸標高')}${this.control('landing.runwayTrackDeg', '真方位（進入）', '°')}</div><div class="glp-two">${this.control('landing.final.pathAngleDeg', '目標パス角', '°')}${this.control('landing.final.groundLengthM', '直線ファイナル長', 'm')}</div><p class="glp-tip">最終パス角は目標値であり、スポイラー・フレアは未モデル化です。</p><div class="glp-preset"><button class="glp-btn glp-btn-soft" data-action="preset-left">左場周プリセット</button><button class="glp-btn glp-btn-soft" data-action="preset-right">右場周プリセット</button></div></div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot blue"></span>離脱側レグ（順方向）</h3><small>${this.project.departureLegs.length}件</small></div>${this.legList('departure')}</div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot green"></span>着陸側レグ（逆方向編集）</h3><small>${this.project.arrivalLegs.length}件</small></div>
   ${this.legList('arrival')}
   <p class="glp-tip">新しいレグは最上流へ追加され、一覧は実際の飛行順です。</p></div>
   <div class="glp-panel"><h3><span class="glp-dot purple"></span>滑らか接続（5次Bezier）</h3>
   ${this.checkbox('display.showHandles', '地図に制御ハンドルを表示')}
   ${this.project.connector ? this.select('connector.adjustmentMode', '接続曲線の調整', [['auto', '自動調整（標準）'], ['manual', '手動調整（ハンドル操作）']]) : `<p class="glp-muted">${this.project.arrivalLegs.length ? '接続がオフになっています。' : '未設計：着陸側レグまたは場周プリセットを追加すると自動生成します。'}</p>`}
   <div class="glp-adds"><button class="glp-btn glp-btn-primary" data-action="optimize" ${!this.project.arrivalLegs.length && !this.project.connector ? 'disabled' : ''}>曲線を自動調整</button><button class="glp-btn glp-btn-sub" data-action="toggle-connector" ${!this.project.arrivalLegs.length && !this.project.connector ? 'disabled' : ''}>${this.project.connector ? '接続を解除' : '接続を作成'}</button></div><p class="glp-tip">接線と曲率を一致させます。大きいバンク角や風で不成立になる場合があります。</p></div>
   <div class="glp-panel"><h3>地図・表示</h3>${this.checkbox('display.showDistanceRings', '着陸点から距離リング')}${this.checkbox('display.showBearingLines', '方位線')}${this.checkbox('display.showWindArrow', '風向矢印')}${this.checkbox('display.showDownwindMargin', 'ダウンウィンド開始の高度余裕を表示')}${this.checkbox('display.showCheckpointAltitudes', '各チェックポイント高度を地図に表示')}${this.opacityControl('display.referenceLineOpacity', '補助線の不透明度')}<div class="glp-two">${this.control('display.distanceRingIntervalM', '距離リング間隔', 'm', '25')}${this.select('display.bearingLineIntervalDeg', '方位線間隔', [['15', '15°'], ['30', '30°'], ['45', '45°'], ['60', '60°'], ['90', '90°'], ['120', '120°'], ['180', '180°']])}</div>${this.control('display.magneticDeclinationDeg', '磁気偏角（東偏＋／西偏−）', '°', '0.1')}<p class="glp-tip">磁気偏角は補助方位線のみへ反映。飛行計算・風向・滑走路方位は真方位基準です。</p></div>`;
        sidebarScroll.scrollTop = preservedScroll;
    }
    renderToolbar() {
        this.root.querySelector('.glp-toolbar').innerHTML = `<div class="glp-brand"><div class="glp-logo">GLP</div><div><strong>Glider Landing Planner</strong><span>飛行経路シミュレーター / テスト版 v0.1.20</span></div></div><div class="glp-tool-group">
  <input type="text" class="glp-project-name" data-path="name" value="${escapeHtml(this.project.name)}" aria-label="プロジェクト名" maxlength="120">
  <button class="glp-btn glp-btn-quiet" data-action="new">新規</button><button class="glp-btn glp-btn-quiet" data-action="undo" ${!this.history.length ? 'disabled' : ''}>元に戻す</button><button class="glp-btn glp-btn-quiet" data-action="redo" ${!this.future.length ? 'disabled' : ''}>やり直す</button><button class="glp-btn glp-btn-quiet" data-action="export" >JSON保存</button><button class="glp-btn glp-btn-quiet" data-action="import">JSON読込</button><button class="glp-btn glp-btn-primary" data-action="print" >PDF印刷</button></div><div class="glp-save-state">● ${escapeHtml(this.status)}</div>`;
    }
    renderResults() {
        const p = this.plan, unit = this.project.display.altitudeUnit;
        const unplanned = p.status === 'unplanned';
        const downwind = p.checkpoints.find(c => c.checkpointKind === 'downwind-start');
        const connected = p.connector?.feasible && downwind?.marginM !== undefined;
        const missing = unplanned ? '—' : '評価不可';
        const { totalDistanceM, totalTimeS, connectionMarginM } = p.metrics;
        const m = `<div class="glp-stat"><small>離脱高度 MSL</small><strong>${fmtAlt(this.project.departure.altitudeMSLM, unit)}</strong></div><div class="glp-stat"><small>推定総距離</small><strong>${totalDistanceM === null ? missing : (totalDistanceM / 1000).toFixed(2) + ' km'}</strong></div><div class="glp-stat"><small>推定飛行時間</small><strong>${totalTimeS === null ? missing : fmtElapsed(totalTimeS)}</strong></div><div class="glp-stat"><small>高度余裕</small><strong>${connectionMarginM === null ? missing : fmtAlt(connectionMarginM, unit)}</strong></div>`;
        this.root.querySelector('.glp-chart').innerHTML = chartSvg(p, this.project);
        this.root.querySelector('.glp-metrics').innerHTML = m;
        this.root.querySelector('.glp-checkpoints').innerHTML = `<div class="glp-overview-top"><strong>チェックポイント高度</strong><small>滑空場基準高度（着陸標高を0とする相対高度）／MSL高度</small></div><div class="glp-tablewrap">${checkpointTable(p, this.project)}</div><p class="glp-tip">高度余裕は接続点での予想・必要高度の差です。着陸地点基準の高さではありません。</p>${downwind && !connected ? '<p class="glp-tip">ダウンウィンド開始の高度余裕は、接続条件が成立していないため評価不可です。</p>' : ''}`;
        const distinct = new Map();
        for (const notice of p.notices)
            distinct.set(JSON.stringify([notice.code, notice.level, notice.legId, notice.message]), notice.message);
        const top = p.notices.filter(i => i.level === 'error');
        const warnings = p.notices.filter(i => i.level === 'warning');
        this.root.querySelector('.glp-diagnostics').innerHTML = `<div class="glp-overview-top"><strong>計算上の診断</strong><span class="glp-status-chip ${top.length ? 'glp-status-bad' : warnings.length ? 'glp-status-warn' : 'glp-status-note'}">${p.status === 'invalid' ? '不成立の条件あり' : unplanned ? '未設計' : warnings.length ? '計算成立（警告あり）' : '計算成立（近似モデル）'}</span></div>
 <div class="glp-diagnostic-bank">Bezier接続最大必要バンク：<strong>${p.connector?.samples.length ? p.connector.maxBankDeg.toFixed(1) + '°' : p.status === 'invalid' ? '評価不可' : '—'}</strong>（設定上限 ${this.project.maxBankDeg.toFixed(1)}°）</div>${this.notice ? `<div class="glp-inline-message">${escapeHtml(this.notice)}</div>` : ''}<div class="glp-alerts">${top.length ? `<p>エラー ${top.length} 件</p>` : ''}${[...distinct.values()].slice(0, 12).map(t => `<div class="glp-alert">${escapeHtml(t)}</div>`).join('') || '<p>計算エラーは報告されていません。</p>'}${distinct.size > 12 ? `<p>ほか ${distinct.size - 12} 件省略</p>` : ''}</div>`;
        this.updateCheckpointSelection();
    }
    updateCheckpointSelection() {
        for (const row of this.root.querySelectorAll('.glp-checkpoints tr[data-leg-id]')) {
            row.classList.toggle('glp-checkpoint-active', !!this.selected && readLegId(row.dataset.legId) === this.selected);
        }
    }
    printReport() {
        const mapEl = this.root.querySelector('.glp-map');
        if (!mapEl)
            throw new Error('地図が見つかりません');
        const diagnostics = this.root.querySelector('.glp-diagnostics');
        const diagnosticsHtml = ['.glp-diagnostic-bank', '.glp-inline-message', '.glp-alerts']
            .map(selector => diagnostics.querySelector(selector)?.outerHTML ?? '').join('');
        openPrintReport({
            project: this.project,
            plan: this.plan,
            mapEl,
            metricsHtml: this.root.querySelector('.glp-metrics').innerHTML,
            diagnosticsHtml,
            title: 'Glider Landing Planner',
            engineVersion: ENGINE_VERSION
        });
    }
    render() {
        this.renderToolbar();
        this.renderSidebar();
        this.renderResults();
        this.map.redraw(this.project, this.plan, this.selected);
        this.renderMobilePanel();
        for (const button of this.root.querySelectorAll('.glp-map-actions button')) {
            button.classList.toggle('glp-active', button.dataset.action === 'pick-' + this.pick);
        }
    }
    /** Only flight inputs can change the geometric connection. Labels and display settings cannot. */
    optimizationSignature(p) {
        const legMotion = (legs) => legs.map(leg => {
            const base = { guidance: leg.guidance, speedOverrideMps: leg.speedOverrideMps };
            if (leg.type === 'straight')
                return { ...base, type: leg.type, lengthM: leg.lengthM, directionDeg: leg.directionDeg };
            return { ...base, type: leg.type, direction: leg.direction, angleDeg: leg.angleDeg, bankDeg: leg.bankDeg, angleInputMode: leg.angleInputMode, targetDirectionDeg: leg.targetDirectionDeg };
        });
        return JSON.stringify({
            speed: p.defaultSpeedMps, bank: p.maxBankDeg, performance: p.performance, wind: p.wind,
            departure: p.departure, landing: p.landing,
            departureLegs: legMotion(p.departureLegs), arrivalLegs: legMotion(p.arrivalLegs)
        });
    }
    recalculate(p, previous) {
        const preliminary = buildPlan(p);
        if (!p.connector || p.connector.adjustmentMode !== 'auto' || !preliminary.connectorReady)
            return preliminary;
        const modeChanged = previous.connector?.adjustmentMode !== 'auto';
        const shapeChanged = this.optimizationSignature(previous) !== this.optimizationSignature(p);
        if (!modeChanged && !shapeChanged)
            return preliminary;
        const optimized = chooseConnector(preliminary.departure.end, preliminary.arrivalStart, p);
        if (!validHandle(optimized.handleStartM) || !validHandle(optimized.handleEndM))
            throw Error('自動生成したBezierハンドル長が不正です。');
        p.connector.handleStartM = optimized.handleStartM;
        p.connector.handleEndM = optimized.handleEndM;
        return buildPlan(p);
    }
    /** Prepare an isolated candidate before changing history or persisted state. */
    prepareProject(project, plan) {
        validateProject(project);
        const result = plan ?? buildPlan(project);
        rememberResolvedTurnAngles(project, result);
        validateProject(project);
        return result;
    }
    /** Common completion path for validated candidates; history is handled by callers. */
    applyProject(project, plan, status) {
        this.project = project;
        this.plan = plan;
        this.notice = '';
        this.status = status;
        this.render();
        this.scheduleSave();
    }
    edit(mutate) {
        // Reject a late startup autosave when any user edit was attempted.
        this.editRevision++;
        const original = this.project;
        const project = structuredClone(original);
        try {
            mutate(project);
            // First arrival-side leg makes a connection meaningful. Returning to zero
            // arrival legs restores the unplanned state; old imported plans stay untouched.
            if (!original.arrivalLegs.length && project.arrivalLegs.length && !project.connector) {
                project.connector = {
                    type: 'bezier_connector', handleStartM: 70, handleEndM: 70,
                    curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto'
                };
            }
            if (original.arrivalLegs.length && !project.arrivalLegs.length)
                project.connector = null;
            validateProject(project);
            if (JSON.stringify(project) === JSON.stringify(original))
                return false;
            // Name and display changes don't affect the flight result. Checkpoint labels do.
            const flightInputsChanged = JSON.stringify({ ...project, name: undefined, display: undefined }) !==
                JSON.stringify({ ...original, name: undefined, display: undefined });
            const plan = this.prepareProject(project, flightInputsChanged ? this.recalculate(project, original) : this.plan);
            this.history.push(structuredClone(original));
            if (this.history.length > 50)
                this.history.shift();
            this.future = [];
            this.applyProject(project, plan, '未保存');
            return true;
        }
        catch (error) {
            this.notice = '変更を適用できません: ' + (error instanceof Error ? error.message : String(error));
            this.render();
            return false;
        }
    }
    scheduleSave() {
        if (this.saveTimer !== undefined)
            window.clearTimeout(this.saveTimer);
        const revision = this.editRevision;
        this.saveTimer = window.setTimeout(() => {
            this.saveTimer = undefined;
            this.status = '保存中…';
            this.renderToolbar();
            // A later edit must not be marked saved by an older asynchronous request.
            void saveSnapshot(structuredClone(this.project)).then(() => {
                if (this.editRevision !== revision)
                    return;
                this.status = '自動保存済';
                this.renderToolbar();
            }).catch(() => {
                if (this.editRevision !== revision)
                    return;
                this.status = '自動保存失敗（JSONを保存してください）';
                this.renderToolbar();
            });
        }, 1000);
    }
    handleChange(input) {
        const path = input.dataset.path;
        if (!path)
            return;
        let value;
        if (input instanceof HTMLInputElement && input.type === 'checkbox') {
            value = input.checked;
        }
        else if (input instanceof HTMLInputElement && (input.type === 'number' || input.type === 'range')) {
            if (input.value.trim() === '') {
                if (path.endsWith('speedOverrideMps') || path.endsWith('directionDeg')) {
                    value = null;
                }
                else {
                    this.notice = '数値項目を空欄にはできません。';
                    this.render();
                    return;
                }
            }
            else {
                value = Number(input.value) / Number(input.dataset.scale ?? 1);
                if (!Number.isFinite(value)) {
                    this.notice = '数値が不正です。';
                    this.render();
                    return;
                }
            }
        }
        else if (input instanceof HTMLSelectElement && path === 'display.bearingLineIntervalDeg') {
            value = Number(input.value);
        }
        else {
            value = input.value;
        }
        if (path.endsWith('.angleInputMode')) {
            const keys = path.split('.');
            const side = keys[0] === 'departureLegs' ? 'departure' : 'arrival';
            const oldLeg = (side === 'departure' ? this.project.departureLegs : this.project.arrivalLegs)[Number(keys[1])];
            if (oldLeg?.type !== 'turn' || (value !== 'angle' && value !== 'target'))
                return;
            const turn = effectiveTurnValues(this.plan, side, oldLeg);
            if (value === 'target' && turn.boundaryDeg === null) {
                this.notice = '旋回境界の方向を評価できないため、方向指定に切り替えられません。';
                this.renderResults();
                return;
            }
            this.edit(project => {
                const selected = (side === 'departure' ? project.departureLegs : project.arrivalLegs)[Number(keys[1])];
                if (selected.type !== 'turn')
                    return;
                selected.angleDeg = turn.angleDeg;
                selected.angleInputMode = value;
                selected.targetDirectionDeg = value === 'target' ? turn.boundaryDeg : null;
            });
            return;
        }
        this.edit(project => {
            const keys = path.split('.');
            let target = project;
            for (let i = 0; i < keys.length - 1; i++)
                target = target[keys[i]];
            const field = keys[keys.length - 1];
            if (path.endsWith('checkpointLabel') && !value)
                delete target[field];
            else
                target[field] = value;
        });
    }
    handleAction(action, data) {
        // Actions that do not change the flight project.
        if (action === 'fit') {
            this.map.fit(this.project, this.plan);
            return;
        }
        if (action === 'pick-departure' || action === 'pick-landing') {
            this.pick = action === 'pick-departure' ? 'departure' : 'landing';
            this.map.setPick(this.pick);
            this.render();
            return;
        }
        if (action === 'import') {
            this.root.querySelector('.glp-import-file').click();
            return;
        }
        if (action === 'export') {
            try {
                downloadProject(this.project);
                this.status = 'JSONを書き出しました';
                this.renderToolbar();
            }
            catch (error) {
                this.notice = String(error);
                this.renderResults();
            }
            return;
        }
        if (action === 'print') {
            try {
                this.printReport();
            }
            catch (error) {
                this.notice = 'PDF印刷を開始できません: ' + (error instanceof Error ? error.message : String(error));
                this.renderResults();
            }
            return;
        }
        // Undo and redo select an already validated project; do not optimize its saved handles.
        if ((action === 'undo' && this.history.length) || (action === 'redo' && this.future.length)) {
            const source = action === 'undo' ? this.history : this.future;
            const dest = action === 'undo' ? this.future : this.history;
            try {
                const project = structuredClone(source.at(-1));
                const plan = this.prepareProject(project);
                this.editRevision++;
                dest.push(structuredClone(this.project));
                source.pop();
                this.applyProject(project, plan, '未保存');
            }
            catch (error) {
                this.notice = '履歴を復元できません: ' + String(error);
                this.renderResults();
            }
            return;
        }
        if (action === 'new') {
            if (window.confirm('現在の計画を破棄して新規計画を作成しますか？（先にJSON保存を推奨）')) {
                const project = createDemo();
                let plan;
                try {
                    plan = this.prepareProject(project);
                }
                catch (error) {
                    this.notice = '新規計画を作成できません: ' + String(error);
                    this.renderResults();
                    return;
                }
                // A new plan starts with the same unoptimized control handles as first launch.
                // Import/restore/undo continue to preserve their saved connector geometry.
                this.editRevision++;
                this.history.push(structuredClone(this.project));
                if (this.history.length > 50)
                    this.history.shift();
                this.future = [];
                this.selected = null;
                this.pick = null;
                this.map.setPick(null);
                this.mobilePanel = null;
                this.applyProject(project, plan, '未保存');
                this.map.showNewProject(project);
            }
            return;
        }
        if (action === 'select-leg') {
            if (data.id)
                this.selectLeg(data.id);
            return;
        }
        if (action === 'close-mobile-panel') {
            this.mobilePanel = null;
            this.renderMobilePanel();
            return;
        }
        if (action === 'edit-departure') {
            this.openMobilePanel({ kind: 'point', point: 'departure' });
            return;
        }
        if (action === 'open-add') {
            if (data.side === 'departure' || data.side === 'arrival')
                this.openMobilePanel({ kind: 'add', side: data.side });
            return;
        }
        // Leg operations. Arrival legs are stored in flight order, so new legs go first.
        if (action === 'add-straight' || action === 'add-turn') {
            if (data.side !== 'departure' && data.side !== 'arrival')
                return;
            const side = data.side;
            const leg = action === 'add-straight' ? makeStraight() : makeTurn();
            const added = this.edit(project => {
                if (side === 'departure')
                    project.departureLegs.push(leg);
                else
                    project.arrivalLegs.unshift(leg);
            });
            if (!added)
                return;
            this.selectLeg(leg.id, true);
            return;
        }
        if (action === 'delete-leg' || action === 'move-up' || action === 'move-down') {
            if (data.side !== 'departure' && data.side !== 'arrival')
                return;
            const side = data.side;
            const id = data.id;
            if (!id)
                return;
            const changed = this.edit(project => {
                const legs = side === 'departure' ? project.departureLegs : project.arrivalLegs;
                const index = legs.findIndex(leg => leg.id === id);
                if (index < 0)
                    return;
                if (action === 'delete-leg') {
                    legs.splice(index, 1);
                }
                else {
                    const destination = index + (action === 'move-up' ? -1 : 1);
                    if (destination < 0 || destination >= legs.length)
                        return;
                    [legs[index], legs[destination]] = [legs[destination], legs[index]];
                }
            });
            if (!changed)
                return;
            // edit() already rendered the moved/remaining legs. Only a deleted
            // selection needs a second redraw to clear its highlight.
            if (action === 'delete-leg' && this.selected === id) {
                this.selected = null;
                if (this.mobilePanel?.kind === 'leg' && this.mobilePanel.id === id)
                    this.mobilePanel = null;
                this.renderSidebar();
                this.map.redraw(this.project, this.plan, null);
                this.updateCheckpointSelection();
                this.renderMobilePanel();
            }
            return;
        }
        if (action === 'preset-left' || action === 'preset-right') {
            if (this.project.arrivalLegs.length && !window.confirm('着陸側レグをプリセットで置き換えますか？'))
                return;
            const changed = this.edit(project => {
                project.arrivalLegs = trafficPreset(action === 'preset-left' ? 'left' : 'right');
                if (!project.connector)
                    project.connector = { type: 'bezier_connector', handleStartM: 70, handleEndM: 70, curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto' };
            });
            if (changed && isCompactLayout() && this.project.arrivalLegs.length) {
                this.selectLeg(this.project.arrivalLegs[0].id, true);
            }
            return;
        }
        // Connector edits retain nonzero endpoint curvature and the existing optimizer.
        if ((action === 'toggle-connector' || action === 'optimize') && !this.project.arrivalLegs.length && !this.project.connector)
            return;
        if (action === 'toggle-connector') {
            this.edit(project => {
                if (project.connector) {
                    project.connector = null;
                    return;
                }
                const span = Math.hypot(this.plan.arrivalStart.e - this.plan.departure.end.e, this.plan.arrivalStart.n - this.plan.departure.end.n);
                const handle = Number.isFinite(span) ? clampHandle(Math.max(20, span * 0.35)) : 70;
                project.connector = {
                    type: 'bezier_connector', handleStartM: handle, handleEndM: handle,
                    curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto'
                };
            });
            return;
        }
        if (action === 'optimize') {
            if (!this.plan.connectorReady)
                return;
            this.edit(project => {
                const result = chooseConnector(this.plan.departure.end, this.plan.arrivalStart, project);
                if (!validHandle(result.handleStartM) || !validHandle(result.handleEndM))
                    throw Error('自動生成したBezierハンドル長が不正です。');
                project.connector = {
                    type: 'bezier_connector', handleStartM: result.handleStartM, handleEndM: result.handleEndM,
                    curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto'
                };
            });
        }
    }
    async doImport(input) {
        const file = input.files?.[0];
        if (!file)
            return;
        try {
            if (file.size > 2000000)
                throw Error('ファイルサイズ上限は2MBです。');
            const rawText = await file.text();
            const project = parseProject(rawText);
            const plan = this.prepareProject(project);
            if (!window.confirm('JSONの内容で現在のプロジェクトを置き換えますか？'))
                return;
            this.editRevision++;
            this.history.push(structuredClone(this.project));
            if (this.history.length > 50)
                this.history.shift();
            this.future = [];
            this.selected = null;
            this.mobilePanel = null;
            this.applyProject(project, plan, 'JSON読込完了');
            this.map.fit(this.project, this.plan);
        }
        catch (error) {
            this.notice = 'JSON読込エラー: ' + (error instanceof Error ? error.message : String(error));
            this.renderResults();
        }
        finally {
            input.value = '';
        }
    }
}
for (const root of document.querySelectorAll('.glp-app-root')) {
    try {
        new App(root);
    }
    catch (e) {
        root.textContent = 'Glider Landing Planner 起動失敗: ' + String(e);
        console.error(e);
    }
}
