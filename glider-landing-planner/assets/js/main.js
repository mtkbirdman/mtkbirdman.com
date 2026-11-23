import { createDemo, makeStraight, makeTurn, trafficPreset, ENGINE_VERSION } from './model/defaults.js';
import { buildPlan } from './planning/plan.js';
import { chooseConnector } from './planning/connector.js';
import { validateProject, parseProject, downloadProject } from './io/project.js';
import { saveSnapshot, latestSnapshot } from './io/autosave.js';
import { PlannerMap } from './ui/map.js';
import { chartSvg, printMapSvg, fmtAlt } from './ui/chart.js';
const escapeHTML = (val) => String(val ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c] ?? c));
const n = (value, d = 2) => Number.isFinite(value) ? Number(value.toFixed(d)) : 0;
const label = (s) => escapeHTML(s);
const options = (items, selected) => items.map(([v, t]) => `<option value="${v}" ${v === selected ? 'selected' : ''}>${t}</option>`).join('');
const engineNote = '教育・机上検討用。実際の飛行可否、地形・障害物、フレア・スポイラー性能は評価しません。';
class App {
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
    constructor(root) {
        this.root = root;
        root.innerHTML = `<div class="glp-shell"><div class="glp-toolbar"></div><div class="glp-app-body"><aside class="glp-sidebar"></aside><section class="glp-main"><div class="glp-mapbar"><div class="glp-map-actions"><button data-action="pick-departure" class="glp-btn glp-btn-quiet">離脱地点を指定</button><button data-action="pick-landing" class="glp-btn glp-btn-quiet">着陸地点を指定</button><button data-action="fit" class="glp-btn glp-btn-quiet">全体表示</button></div><span class="glp-map-hint">地点・直線端点・紫のハンドルはドラッグできます</span></div><div class="glp-map"></div><div class="glp-overview"><div class="glp-overview-top"><strong>高度プロファイル</strong><span class="glp-chart-legend"><i class="glp-line-blue"></i>予想高度 <i class="glp-line-orange"></i>必要高度</span></div><div class="glp-chart"></div><div class="glp-metrics"></div><div class="glp-checkpoints"></div><div class="glp-diagnostics"></div></div></section></div><div class="glp-footnote">${engineNote}　｜　地図：<a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>（地理院タイル）</div><div class="glp-print"></div><input class="glp-import-file" type="file" accept=".json,.glp.json,application/json" hidden></div>`;
        this.map = new PlannerMap(root.querySelector('.glp-map'), {
            onPick: (type, lat, lon) => { this.pick = null; this.map.setPick(null); this.edit(p => { if (type === 'departure') {
                p.departure.latDeg = lat;
                p.departure.lonDeg = lon;
            }
            else {
                p.landing.point.latDeg = lat;
                p.landing.point.lonDeg = lon;
            } }); },
            onLegLength: (side, id, lengthM) => this.edit(p => { const l = (side === 'departure' ? p.departureLegs : p.arrivalLegs).find(v => v.id === id); if (l?.type === 'straight')
                l.lengthM = lengthM; }),
            onFinalLength: lengthM => this.edit(p => { p.landing.final.groundLengthM = lengthM; }),
            onHandle: (end, lengthM) => this.edit(p => { if (p.connector) {
                if (end === 'start')
                    p.connector.handleStartM = lengthM;
                else
                    p.connector.handleEndM = lengthM;
            } })
        });
        root.addEventListener('click', e => { const btn = e.target.closest('button[data-action]'); if (btn) {
            e.preventDefault();
            this.handleAction(btn.dataset.action ?? '', btn.dataset);
        } });
        root.addEventListener('change', e => {
            const el = e.target;
            if (el.matches('[data-path]'))
                this.handleChange(el);
            else if (el.classList.contains('glp-import-file'))
                void this.doImport(el);
        });
        this.render();
        this.map.fit(this.project, this.plan);
        void latestSnapshot().then(p => {
            if (!p)
                return;
            const valid = validateProject(p);
            if (window.confirm('前回の自動保存プロジェクトがあります。復元しますか？')) {
                this.project = valid;
                this.plan = buildPlan(valid);
                this.history = [];
                this.future = [];
                this.status = '自動保存から復元';
                this.render();
                this.map.fit(this.project, this.plan);
            }
        }).catch(() => { });
    }
    getField(path) { let a = this.project; for (const k of path.split('.'))
        a = a[k]; return a; }
    control(path, title, unit = '', step = 'any', scale = 1) {
        const raw = this.getField(path);
        const value = typeof raw === 'number' ? n(raw * scale, scale === 3.6 ? 2 : 5) : raw ?? '';
        return `<label class="glp-field"><span>${title}${unit ? ` <em>${unit}</em>` : ''}</span><input type="number" step="${step}" value="${escapeHTML(value)}" data-path="${path}" data-scale="${scale}"></label>`;
    }
    checkbox(path, title) { return `<label class="glp-check"><input type="checkbox" data-path="${path}" ${this.getField(path) ? 'checked' : ''}><span>${title}</span></label>`; }
    select(path, title, items) { return `<label class="glp-field"><span>${title}</span><select data-path="${path}">${options(items, String(this.getField(path)))}</select></label>`; }
    legList(side) {
        const legs = side === 'departure' ? this.project.departureLegs : this.project.arrivalLegs;
        const modeOpts = [['heading', '機首方位指定（風で偏流）'], ['ground-track', '対地航跡指定（風を補正）']];
        const lines = legs.map((leg, i) => {
            const active = this.selected === leg.id;
            const title = leg.type === 'straight' ? `直線 ${Math.round(leg.lengthM)} m` : `${leg.direction === 'right' ? '右' : '左'}旋回 ${leg.angleDeg}° / バンク${leg.bankDeg}°`;
            const group = side === 'departure' ? 'departureLegs' : 'arrivalLegs';
            const base = `${group}.${i}`;
            return `<div class="glp-leg ${active ? 'glp-leg-active' : ''}"><div class="glp-leg-line"><button class="glp-leg-title" data-action="select-leg" data-id="${leg.id}" title="編集"><span class="glp-leg-symbol">${leg.type === 'straight' ? '━' : '⤷'}</span><span>${i + 1}. ${label(title)}</span></button><span class="glp-leg-buttons"><button title="上に移動" data-action="move-up" data-side="${side}" data-id="${leg.id}">↑</button><button title="下に移動" data-action="move-down" data-side="${side}" data-id="${leg.id}">↓</button><button title="削除" data-action="delete-leg" data-side="${side}" data-id="${leg.id}">×</button></span></div>
     ${active ? `<div class="glp-leg-edit"><div class="glp-two">${this.select(base + '.guidance', '風補正方式', modeOpts)}${this.control(base + '.speedOverrideMps', '個別速度（空欄＝共通）', 'km/h', '0.1', 3.6)}</div>
       ${leg.type === 'straight' ? `<div class="glp-two">${this.control(base + '.lengthM', '直線距離', 'm')}${this.control(base + '.directionDeg', '指定方位（空欄＝継承）', '°')}</div>` : `<div class="glp-two">${this.select(base + '.direction', '旋回方向', [['right', '右'], ['left', '左']])}${this.control(base + '.angleDeg', '旋回角', '°')}</div><div class="glp-two">${this.control(base + '.bankDeg', 'バンク角', '°')}</div>`}
       <label class="glp-field"><span>チェックポイント名</span><input type="text" maxlength="70" value="${escapeHTML(leg.checkpointLabel ?? '')}" data-path="${base}.checkpointLabel"></label></div>` : ''}</div>`;
        }).join('');
        return `<div class="glp-leg-list">${lines || '<p class="glp-muted">まだレグがありません。</p>'}</div><div class="glp-adds"><button class="glp-btn glp-btn-sub" data-action="add-straight" data-side="${side}">＋ 直線</button><button class="glp-btn glp-btn-sub" data-action="add-turn" data-side="${side}">＋ 旋回</button></div>`;
    }
    renderSidebar() {
        this.root.querySelector('.glp-sidebar').innerHTML = `
   <div class="glp-panel-heading"><div class="glp-eyebrow">FLIGHT PARAMETERS</div><h2>飛行条件・経路設計</h2></div>
   <div class="glp-panel"><h3>基本条件</h3><div class="glp-two">${this.control('defaultSpeedMps', '共通対気速度', 'km/h', '0.1', 3.6)}${this.control('maxBankDeg', '最大バンク', '°')}</div>
   ${this.select('performance.model', '沈下計算モデル', [['fixed-ld', 'A：固定L/D'], ['quadratic-drag', 'B：簡易抗力モデル']])}
   <div class="glp-two">${this.control('performance.fixedLD', 'モデルA：L/D')}${this.control('performance.bestLD', 'モデルB：最大L/D')}</div>${this.control('performance.bestLDAtMps', 'モデルB：最大L/D速度', 'km/h', '0.1', 3.6)}
   <div class="glp-two">${this.control('wind.fromDeg', '風向（吹いてくる）', '°')}${this.control('wind.speedMps', '風速', 'm/s')}</div></div>
   <div class="glp-panel"><h3><span class="glp-dot blue"></span>離脱地点・条件</h3><div class="glp-two">${this.control('departure.latDeg', '緯度', '°')}${this.control('departure.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('departure.altitudeMSLM', '離脱高度 MSL', 'm')}${this.control('departure.headingDeg', '開始機首方位', '°')}</div><p class="glp-tip">地図の「離脱地点を指定」または青いマーカーのドラッグで変更できます。</p></div>
   <div class="glp-panel"><h3><span class="glp-dot green"></span>着陸地点・ファイナル</h3><div class="glp-two">${this.control('landing.point.latDeg', '緯度', '°')}${this.control('landing.point.lonDeg', '経度', '°')}</div><div class="glp-two">${this.control('landing.elevationMSLM', '着陸標高 MSL', 'm')}${this.control('landing.runwayTrackDeg', '真方位（進入）', '°')}</div><div class="glp-two">${this.control('landing.final.pathAngleDeg', '目標パス角', '°')}${this.control('landing.final.groundLengthM', '直線ファイナル長', 'm')}</div><p class="glp-tip">最終パス角は目標値であり、スポイラー・フレアは未モデル化です。</p></div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot blue"></span>離脱側レグ（順方向）</h3><small>${this.project.departureLegs.length}件</small></div>${this.legList('departure')}</div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot green"></span>着陸側レグ（逆方向編集）</h3><small>${this.project.arrivalLegs.length}件</small></div>
   <div class="glp-preset"><button class="glp-btn glp-btn-soft" data-action="preset-left">左場周プリセット</button><button class="glp-btn glp-btn-soft" data-action="preset-right">右場周プリセット</button></div>${this.legList('arrival')}
   <p class="glp-tip">新しいレグは最上流へ追加され、一覧は実際の飛行順です。</p></div>
   <div class="glp-panel"><h3><span class="glp-dot purple"></span>滑らか接続（5次Bezier）</h3>
   ${this.checkbox('display.showHandles', '地図に制御ハンドルを表示')}
   ${this.project.connector ? `<div class="glp-two">${this.control('connector.handleStartM', '始端ハンドル', 'm')}${this.control('connector.handleEndM', '終端ハンドル', 'm')}</div>` : '<p class="glp-muted">接続がオフになっています。</p>'}
   <div class="glp-adds"><button class="glp-btn glp-btn-primary" data-action="optimize">曲線を自動調整</button><button class="glp-btn glp-btn-sub" data-action="toggle-connector">${this.project.connector ? '接続を解除' : '接続を作成'}</button></div><p class="glp-tip">接線と曲率を一致させます。大きいバンク角や風で不成立になる場合があります。</p></div>
   <div class="glp-panel"><h3>地図・表示</h3>${this.select('display.altitudeUnit', '高度表示単位', [['ft', 'ft'], ['m', 'm']])}${this.checkbox('display.showDistanceRings', '着陸点から距離リング')}${this.checkbox('display.showBearingLines', '方位線')}${this.checkbox('display.showWindArrow', '風向矢印')}</div>`;
    }
    renderToolbar() {
        this.root.querySelector('.glp-toolbar').innerHTML = `<div class="glp-brand"><div class="glp-logo">GLP</div><div><strong>Glider Landing Planner</strong><span>飛行経路シミュレーター / テスト版 v0.1.1</span></div></div><div class="glp-tool-group">
  <input type="text" class="glp-project-name" data-path="name" value="${escapeHTML(this.project.name)}" aria-label="プロジェクト名" maxlength="120">
  <button class="glp-btn glp-btn-quiet" data-action="new">新規</button><button class="glp-btn glp-btn-quiet" data-action="undo" ${!this.history.length ? 'disabled' : ''}>元に戻す</button><button class="glp-btn glp-btn-quiet" data-action="redo" ${!this.future.length ? 'disabled' : ''}>やり直す</button><button class="glp-btn glp-btn-quiet" data-action="export">JSON保存</button><button class="glp-btn glp-btn-quiet" data-action="import">JSON読込</button><button class="glp-btn glp-btn-primary" data-action="print">PDF印刷</button></div><div class="glp-save-state">● ${escapeHTML(this.status)}</div>`;
    }
    renderResults() {
        const p = this.plan, unit = this.project.display.altitudeUnit;
        const gap = p.checkpoints.find(c => c.label === '自動接続 終点') ?? p.checkpoints[0];
        const m = `<div class="glp-stat"><small>離脱高度 MSL</small><strong>${fmtAlt(this.project.departure.altitudeMSLM, unit)}</strong></div><div class="glp-stat"><small>推定総距離</small><strong>${(p.forwardAll.at(-1)?.sM ?? 0) / 1000 < 100 ? ((p.forwardAll.at(-1)?.sM ?? 0) / 1000).toFixed(2) : '—'} km</strong></div><div class="glp-stat"><small>接続区間の高度差</small><strong>${gap ? fmtAlt(gap.marginM, unit) : '—'}</strong></div><div class="glp-stat"><small>接続最大必要バンク</small><strong>${p.connector ? p.connector.maxBankDeg.toFixed(1) + '°' : '—'}</strong></div>`;
        this.root.querySelector('.glp-chart').innerHTML = chartSvg(p, this.project);
        this.root.querySelector('.glp-metrics').innerHTML = m;
        this.root.querySelector('.glp-checkpoints').innerHTML = `<div class="glp-overview-top"><strong>チェックポイント高度余裕</strong><small>予想到達高度 − 逆算必要高度</small></div><div class="glp-tablewrap"><table class="glp-table"><thead><tr><th>地点</th><th>予想</th><th>必要</th><th>差</th></tr></thead><tbody>${p.checkpoints.map(c => `<tr><td>${label(c.label)}</td><td>${fmtAlt(c.predictedM, unit)}</td><td>${fmtAlt(c.requiredM, unit)}</td><td class="${c.marginM < 0 ? 'glp-negative' : 'glp-positive'}">${c.marginM > 0 ? '+' : ''}${fmtAlt(c.marginM, unit)}</td></tr>`).join('')}</tbody></table></div>`;
        const distinct = new Map();
        for (const notice of p.notices)
            distinct.set(notice.code + notice.message, notice.message);
        const top = p.notices.filter(i => i.level === 'error');
        const warnings = p.notices.filter(i => i.level === 'warning');
        this.root.querySelector('.glp-diagnostics').innerHTML = `<div class="glp-overview-top"><strong>計算上の診断</strong><span class="glp-status-chip ${top.length ? 'glp-status-bad' : warnings.length ? 'glp-status-warn' : 'glp-status-note'}">${top.length ? '不成立の条件あり' : warnings.length ? '警告あり' : '近似モデル'}</span></div>
 ${this.notice ? `<div class="glp-inline-message">${label(this.notice)}</div>` : ''}<div class="glp-alerts">${[...distinct.values()].slice(0, 12).map(t => `<div class="glp-alert">${label(t)}</div>`).join('') || '<p>計算エラーは報告されていません。</p>'}</div>`;
    }
    renderPrint() {
        const p = this.plan, pr = this.project, unit = pr.display.altitudeUnit;
        const container = this.root.querySelector('.glp-print');
        container.innerHTML = `<div class="glp-print-head"><h1>Glider Landing Planner ｜ 飛行経路計画</h1><p>${label(pr.name)} ｜ エンジン ${ENGINE_VERSION} ｜ 出力日 ${new Date().toLocaleString('ja-JP')}</p></div><div class="glp-print-figure">${printMapSvg(p)}</div><div class="glp-print-graph">${chartSvg(p, pr)}</div>
 <div class="glp-print-data"><div>離脱高度：${fmtAlt(pr.departure.altitudeMSLM, unit)}</div><div>着陸標高：${fmtAlt(pr.landing.elevationMSLM, unit)}</div><div>速度：${(pr.defaultSpeedMps * 3.6).toFixed(1)} km/h</div><div>風：${pr.wind.fromDeg}° / ${pr.wind.speedMps} m/s</div><div>モデル：${pr.performance.model === 'fixed-ld' ? 'A 固定L/D' : 'B 簡易抗力'}</div><div>目標進入角：${pr.landing.final.pathAngleDeg}°</div></div>
 <table class="glp-print-table"><thead><tr><th>チェックポイント</th><th>予想</th><th>必要</th><th>差</th></tr></thead><tbody>${p.checkpoints.map(c => `<tr><td>${label(c.label)}</td><td>${fmtAlt(c.predictedM, unit)}</td><td>${fmtAlt(c.requiredM, unit)}</td><td>${fmtAlt(c.marginM, unit)}</td></tr>`).join('')}</tbody></table><p class="glp-print-warning">${engineNote}　計算の高度余裕は安全余裕を保証しません。高度はMSL基準。モデルBも実機ポーラーではありません。</p><p class="glp-print-source">経路図はローカルENU座標から描画（背景地図なし）。オンライン地図表示には国土地理院「地理院タイル（淡色地図）」を使用しています。<br>警告：${p.notices.filter(n => n.level !== 'info').slice(0, 5).map(n => label(n.message)).join(' / ') || '特記事項なし'}</p>`;
    }
    render() {
        this.renderToolbar();
        this.renderSidebar();
        this.renderResults();
        this.renderPrint();
        this.map.redraw(this.project, this.plan);
        const bb = this.root.querySelectorAll('.glp-map-actions button');
        bb.forEach(b => b.classList.toggle('glp-active', b.dataset.action === 'pick-' + this.pick));
    }
    edit(mutate) {
        const p = structuredClone(this.project), original = structuredClone(this.project);
        try {
            mutate(p);
            validateProject(p);
            const next = buildPlan(p);
            this.history.push(original);
            if (this.history.length > 50)
                this.history.shift();
            this.future = [];
            this.project = p;
            this.plan = next;
            this.notice = '';
            this.status = '未保存';
            this.render();
            this.scheduleSave();
        }
        catch (e) {
            this.notice = '変更を適用できません: ' + (e instanceof Error ? e.message : String(e));
            this.renderResults();
        }
    }
    scheduleSave() { if (this.saveTimer)
        window.clearTimeout(this.saveTimer); this.saveTimer = window.setTimeout(() => { this.status = '保存中…'; this.renderToolbar(); void saveSnapshot(this.project).then(() => { this.status = '自動保存済'; this.renderToolbar(); }).catch(() => { this.status = '自動保存失敗（JSONを保存してください）'; this.renderToolbar(); }); }, 1000); }
    handleChange(input) {
        const path = input.dataset.path;
        if (!path)
            return;
        let value;
        if (input instanceof HTMLInputElement && input.type === 'checkbox')
            value = input.checked;
        else if (input instanceof HTMLInputElement && input.type === 'number') {
            if (input.value.trim() === '') {
                if (path.endsWith('speedOverrideMps') || path.endsWith('directionDeg'))
                    value = null;
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
        else
            value = input.value;
        this.edit(p => { const keys = path.split('.'); let target = p; for (let i = 0; i < keys.length - 1; i++)
            target = target[keys[i]]; if (path.endsWith('checkpointLabel') && !value)
            delete target[keys.at(-1)];
        else
            target[keys.at(-1)] = value; });
    }
    handleAction(action, data) {
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
            catch (e) {
                this.notice = String(e);
                this.renderResults();
            }
            return;
        }
        if (action === 'print') {
            this.renderPrint();
            window.print();
            return;
        }
        if (action === 'undo') {
            if (this.history.length) {
                this.future.push(structuredClone(this.project));
                this.project = this.history.pop();
                this.plan = buildPlan(this.project);
                this.render();
                this.scheduleSave();
            }
            return;
        }
        if (action === 'redo') {
            if (this.future.length) {
                this.history.push(structuredClone(this.project));
                this.project = this.future.pop();
                this.plan = buildPlan(this.project);
                this.render();
                this.scheduleSave();
            }
            return;
        }
        if (action === 'new') {
            if (window.confirm('現在の計画を破棄して新規計画を作成しますか？（先にJSON保存を推奨）')) {
                this.edit(p => Object.assign(p, createDemo()));
                this.map.fit(this.project, this.plan);
            }
            return;
        }
        if (action === 'select-leg') {
            this.selected = data.id ?? null;
            this.renderSidebar();
            return;
        }
        if (action === 'add-straight' || action === 'add-turn') {
            const leg = action === 'add-straight' ? makeStraight() : makeTurn();
            const side = data.side;
            this.edit(p => { if (side === 'departure')
                p.departureLegs.push(leg);
            else
                p.arrivalLegs.unshift(leg); });
            this.selected = leg.id;
            this.renderSidebar();
            return;
        }
        if (action === 'delete-leg' || action === 'move-up' || action === 'move-down') {
            const side = data.side, id = data.id;
            this.edit(p => {
                const list = side === 'departure' ? p.departureLegs : p.arrivalLegs;
                const i = list.findIndex(l => l.id === id);
                if (i < 0)
                    return;
                if (action === 'delete-leg')
                    list.splice(i, 1);
                else {
                    const dest = i + (action === 'move-up' ? -1 : 1);
                    if (dest < 0 || dest >= list.length)
                        return;
                    [list[i], list[dest]] = [list[dest], list[i]];
                }
            });
            if (action === 'delete-leg' && this.selected === id)
                this.selected = null;
            this.renderSidebar();
            return;
        }
        if (action === 'preset-left' || action === 'preset-right') {
            if (this.project.arrivalLegs.length && !window.confirm('着陸側レグをプリセットで置き換えますか？'))
                return;
            this.edit(p => { p.arrivalLegs = trafficPreset(action === 'preset-left' ? 'left' : 'right'); });
            return;
        }
        if (action === 'toggle-connector') {
            this.edit(p => { p.connector = p.connector ? null : { type: 'bezier_connector', handleStartM: Math.max(20, Math.hypot(this.plan.arrivalStart.e - this.plan.departure.end.e, this.plan.arrivalStart.n - this.plan.departure.end.n) * .35), handleEndM: Math.max(20, Math.hypot(this.plan.arrivalStart.e - this.plan.departure.end.e, this.plan.arrivalStart.n - this.plan.departure.end.n) * .35), curvaturePolicy: 'match-neighbors' }; });
            return;
        }
        if (action === 'optimize') {
            const result = chooseConnector(this.plan.departure.end, this.plan.arrivalStart, this.project);
            this.edit(p => { p.connector = { type: 'bezier_connector', handleStartM: result.handleStartM, handleEndM: result.handleEndM, curvaturePolicy: 'match-neighbors' }; });
            return;
        }
    }
    async doImport(input) {
        const file = input.files?.[0];
        if (!file)
            return;
        try {
            if (file.size > 2_000_000)
                throw Error('ファイルサイズ上限は2MBです。');
            const p = parseProject(await file.text());
            if (!window.confirm('JSONの内容で現在のプロジェクトを置き換えますか？'))
                return;
            this.history.push(structuredClone(this.project));
            this.project = p;
            this.plan = buildPlan(p);
            this.future = [];
            this.selected = null;
            this.status = 'JSON読込完了';
            this.render();
            this.map.fit(this.project, this.plan);
            this.scheduleSave();
        }
        catch (e) {
            this.notice = 'JSON読込エラー: ' + (e instanceof Error ? e.message : String(e));
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
