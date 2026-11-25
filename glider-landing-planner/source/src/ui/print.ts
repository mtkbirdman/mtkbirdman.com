import type { PlanResult, ProjectV1, FlightLeg } from '../model/types.js';
import { chartSvg, fmtAlt } from './chart.js';
import { checkpointTable } from './checkpoints.js';
import { escapeHtml } from './escape-html.js';
interface PrintRequest {
    project: ProjectV1;
    plan: PlanResult;
    mapEl: HTMLElement;
    metricsHtml: string;
    diagnosticsHtml: string;
    title: string;
    engineVersion: string;
}
const num = (value: number, digits = 1) => Number(value.toFixed(digits)).toLocaleString('ja-JP');
const speed = (value: number) => `${num(value * 3.6)} km/h`;
const legLabel = (leg: FlightLeg) => leg.type === 'straight' ? `直線 ${num(leg.lengthM, 0)} m` : `${leg.direction === 'right' ? '右' : '左'}旋回 ${num(leg.angleDeg, 0)}°、バンク ${num(leg.bankDeg, 0)}°`;
/**
 * Uses the ACTUAL rendered Leaflet panes (GSI tile IMG elements, SVG routes,
 * labels, rings and bearings). Cloning avoids CORS canvas taint and preserves
 * the same viewport as the editor. A new document avoids Cocoon print layout
 * and the empty pages caused by visibility:hidden on WordPress ancestors.
 */
export function printableMapHtml(mapEl: HTMLElement): {
    html: string;
    width: number;
    height: number;
    tileCount: number;
} {
    const width = Math.max(1, mapEl.clientWidth || mapEl.offsetWidth || 960);
    const height = Math.max(1, mapEl.clientHeight || mapEl.offsetHeight || 510);
    const cloned = mapEl.cloneNode(true) as HTMLElement;
    // Map editor buttons don't belong in the report, but keep attribution.
    for (const el of cloned.querySelectorAll('.leaflet-control-zoom,.leaflet-control-scale'))
        el.remove();
    const tileCount = cloned.querySelectorAll('.leaflet-tile-pane img.leaflet-tile').length;
    // Use the exact pixel geometry from the live map; scale the whole copy uniformly.
    // 995px corresponds to approximately 263mm at CSS 96dpi (A4 landscape with margins).
    const scale = Math.min(1, 995 / width);
    const renderedHeight = Math.ceil(height * scale);
    cloned.removeAttribute('id');
    cloned.classList.add('glp-report-live-map');
    cloned.style.width = `${width}px`;
    cloned.style.height = `${height}px`;
    cloned.style.transform = `scale(${scale})`;
    cloned.style.transformOrigin = 'top left';
    cloned.style.position = 'relative';
    return {
        html: `<div class="glp-report-map-frame" style="width:${Math.ceil(width * scale)}px;max-width:100%;height:${renderedHeight}px"><div class="glp-report-map-viewport">${cloned.outerHTML}</div></div>`,
        width, height, tileCount
    };
}
function legTable(title: string, legs: FlightLeg[], defaultSpeedMps: number, solvedAngles: Map<string, number>): string {
    if (!legs.length)
        return '';
    const rows = legs.map((leg, i) => `<tr><td>${i + 1}</td><td>${escapeHtml(legLabel(leg.type === 'turn' ? { ...leg, angleDeg: solvedAngles.get(leg.id) ?? leg.angleDeg } : leg))}</td><td>${escapeHtml(leg.guidance === 'heading' ? '機首方位指定' : '対地航跡指定')}</td><td>${speed(leg.speedOverrideMps ?? defaultSpeedMps)}</td></tr>`).join('');
    return `<section class="glp-report-section glp-report-section-legs"><h2>${escapeHtml(title)}</h2><table class="glp-report-table"><thead><tr><th>No.</th><th>内容</th><th>飛行方式</th><th>速度</th></tr></thead><tbody>${rows}</tbody></table></section>`;
}
export function reportHtml(req: PrintRequest): string {
    const { project: p, plan, metricsHtml, diagnosticsHtml, title, engineVersion } = req;
    const map = printableMapHtml(req.mapEl);
    const unit = p.display.altitudeUnit;
    const solvedAngles = new Map([...plan.departure.segments, ...plan.arrival.segments].filter(s => s.leg.type === 'turn').map(s => [s.leg.id, s.leg.type === 'turn' ? s.leg.angleDeg : 0]));
    const attrs = [
        ['離脱高度（MSL）', fmtAlt(p.departure.altitudeMSLM, unit)],
        ['着陸標高（MSL）', fmtAlt(p.landing.elevationMSLM, unit)],
        ['基準対気速度', speed(p.defaultSpeedMps)],
        ['最大バンク', `${p.maxBankDeg}°`],
        ['風向・風速', `${p.wind.fromDeg}°から / ${p.wind.speedMps} m/s`],
        ['沈下モデル', p.performance.model === 'quadratic-drag' ? `B：簡易抗力（最大L/D ${p.performance.quadratic.bestLD}、最良滑空速度 ${speed(p.performance.quadratic.bestLDAtMps)}）` : `C：速度別ポーラー（${p.performance.polarTable?.points.length ?? 0}点、${p.performance.polarTable ? speed(p.performance.polarTable.points[0].speedMps) + '〜' + speed(p.performance.polarTable.points.at(-1)!.speedMps) : '未設定'}）`],
        ['目標ファイナルパス角', `${p.landing.final.pathAngleDeg}°`],
        ['直線ファイナル長', `${num(p.landing.final.groundLengthM, 0)} m`],
        ['滑走路真方位', `${p.landing.runwayTrackDeg}°`],
        ['磁気偏角（東偏＋）', `${p.display.magneticDeclinationDeg}°`],
        ['高度表示単位', unit]
    ];
    const conditions = attrs.map(([key, value]) => `<div class="glp-report-data-item"><small>${escapeHtml(key)}</small><strong>${escapeHtml(value)}</strong></div>`).join('');
    const leafletCss = (document.querySelector<HTMLLinkElement>('link[href*="leaflet.css"]')?.href ?? 'https://unpkg.com/leaflet@1.9.4/dist/leaflet.css');
    const printCss = new URL('../../css/print.css', import.meta.url).href + '?ver=0.1.16';
    const time = new Date().toLocaleString('ja-JP');
    return `<!doctype html><html lang="ja"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>${escapeHtml(p.name)} - ${escapeHtml(title)} PDF印刷</title><link rel="stylesheet" href="${escapeHtml(leafletCss)}"><link rel="stylesheet" href="${escapeHtml(printCss)}"></head>
 <body class="glp-report"><div class="glp-report-screen-toolbar"><span>印刷専用プレビュー｜地図のタイルが読み込まれたら「印刷 / PDF保存」を選択してください。</span><button type="button" onclick="window.print()">印刷 / PDF保存</button></div>
 <main class="glp-report-document">
 <header class="glp-report-head"><div><h1>${escapeHtml(title)}</h1><p>飛行経路計画（教育・机上検討用）</p></div><p>${escapeHtml(p.name)}<br>エンジン ${escapeHtml(engineVersion)}｜${escapeHtml(time)}</p></header>
 <section class="glp-report-section glp-report-map-section"><h2>飛行経路図</h2>${map.html}<p class="glp-report-caption">地図：国土地理院「地理院タイル（淡色地図）」に飛行経路・高度情報等を追記。<br>地図は現在の画面表示範囲を使用しています。印刷前に地図の縮尺・範囲を調整してください。</p></section>
 <section class="glp-report-section glp-report-chart-section"><h2>高度プロファイル</h2><div class="glp-report-chart">${chartSvg(plan, p)}</div><p class="glp-report-caption">青：離脱側予想高度（Bezier接続を含む）。緑：着陸側の逆算必要高度。横軸は推定経過飛行時間。</p></section>
 <section class="glp-report-section"><h2>計算結果</h2><div class="glp-report-metrics">${metricsHtml}</div></section>
 <section class="glp-report-section"><h2>飛行条件</h2><div class="glp-report-data">${conditions}</div></section>
 <section class="glp-report-section"><h2>チェックポイント高度</h2>${checkpointTable(plan, p, true)}</section>
 ${legTable('離脱側レグ（飛行順）', p.departureLegs, p.defaultSpeedMps, solvedAngles)}
 ${legTable('着陸側レグ（飛行順・逆方向編集）', p.arrivalLegs, p.defaultSpeedMps, solvedAngles)}
 <section class="glp-report-section"><h2>計算上の診断</h2><div class="glp-report-diagnostics">${diagnosticsHtml}</div></section>
 <footer class="glp-report-footer">教育・机上検討用。実際の飛行可否・地形・障害物・フレア・スポイラー性能は評価していません。計算上の高度余裕は安全余裕を保証しません。<br>地図出典：国土地理院（地理院タイル） https://maps.gsi.go.jp/development/ichiran.html</footer>
 </main></body></html>`;
}
export function openPrintReport(req: PrintRequest): void {
    // Must run synchronously during the button click to avoid popup blockers.
    const popup = window.open('', '_blank');
    if (!popup)
        throw new Error('印刷プレビューがブロックされました。ブラウザのポップアップ許可を確認してください。');
    const report = reportHtml(req);
    popup.document.open();
    popup.document.write(report);
    popup.document.close();
    popup.focus();
    // Wait for tile images and styles before invoking the print dialog.
    // If a tile failed, keep the preview open for inspection and manual retry.
    void (async () => {
        try {
            // Printing before CSS arrives can hide Leaflet tiles and vector overlays.
            const resources = [
                ...Array.from(popup.document.images).map(img => img.complete ? Promise.resolve() : new Promise<void>(resolve => {
                    img.addEventListener('load', () => resolve(), { once: true });
                    img.addEventListener('error', () => resolve(), { once: true });
                })),
                ...Array.from(popup.document.querySelectorAll<HTMLLinkElement>('link[rel="stylesheet"]')).map(link => link.sheet ? Promise.resolve() : new Promise<void>(resolve => {
                    link.addEventListener('load', () => resolve(), { once: true });
                    link.addEventListener('error', () => resolve(), { once: true });
                }))
            ];
            await Promise.race([Promise.all(resources), new Promise<void>(resolve => setTimeout(resolve, 7500))]);
            if (popup.closed)
                return;
            const visibleTiles = Array.from(popup.document.querySelectorAll<HTMLImageElement>('.leaflet-tile-pane img.leaflet-tile'));
            const loadedTiles = visibleTiles.filter(img => img.naturalWidth > 0);
            if (!loadedTiles.length) {
                const hint = popup.document.querySelector('.glp-report-screen-toolbar span');
                if (hint)
                    hint.textContent = '地図タイルが読み込まれていません。ネットワークを確認し、読み込み後に印刷してください。';
                return; // Do not automatically print a blank basemap.
            }
            await new Promise<void>(resolve => setTimeout(resolve, 450));
            popup.print();
        }
        catch { /* the visible print button is a manual fallback */ }
    })();
}
