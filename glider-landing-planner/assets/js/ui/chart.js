import { heightForDisplay } from './height-units.js?v=0.1.21';
export const fmtAlt = (m, u) => `${Math.round(heightForDisplay(m, u)).toLocaleString('ja-JP')} ${u}`;
export const fmtElapsed = (seconds) => {
    const s = Math.max(0, Math.round(seconds));
    return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
};
/**
 * Use genuine elapsed seconds, not ENU distance or a distance-based display offset.
 * Arrival samples already carry the canonical time offset assigned by buildPlan.
 * Do not recompute a second time axis here.
 * The two curves intentionally remain separate even when their x coordinates meet.
 */
export function altitudeTimeSeries(plan) {
    const prefix = (samples) => { const i = samples.findIndex(s => s.altitudeValid === false || !Number.isFinite(s.tS) || !Number.isFinite(s.h)); return i < 0 ? samples : samples.slice(0, i); };
    const departure = prefix([...plan.departure.samples, ...(plan.connector?.feasible ? plan.connector.samples.slice(1) : [])]);
    // Reverse integration propagates invalid altitude towards the beginning of the arrival.
    const raw = plan.requiredAll;
    let first = 0;
    while (first < raw.length && (raw[first].altitudeValid === false || !Number.isFinite(raw[first].tS)))
        first++;
    const arrival = prefix(raw.slice(first));
    const durationS = plan.metrics.totalTimeS ?? Math.max(0, departure.at(-1)?.tS ?? 0, arrival.at(-1)?.tS ?? 0);
    return { departure, arrival, durationS, timingProvisional: plan.status !== 'valid' };
}
export function chartSvg(plan, project, width = 780, height = 200) {
    // An unplanned new project has no continuous route or meaningful flight-time axis.
    if (plan.status === 'unplanned')
        return `<svg viewBox="0 0 ${width} ${height}" role="img" aria-label="経路は未設計です"><text x="${width / 2}" y="${height / 2}" text-anchor="middle" fill="#64748b" font-size="15">着陸側レグまたは場周プリセットの追加後に高度プロファイルを表示します</text></svg>`;
    // Blue: forecast including quintic Bezier connector. Green: independent arrival requirement.
    const { departure, arrival, durationS, timingProvisional } = altitudeTimeSeries(plan);
    const all = [...departure, ...arrival];
    if (!all.length)
        return '<svg viewBox="0 0 780 200"></svg>';
    const left = 52, right = 16, top = 23, bottom = 35;
    let ymin = Infinity, ymax = -Infinity;
    for (const v of all) {
        if (v.h < ymin)
            ymin = v.h;
        if (v.h > ymax)
            ymax = v.h;
    }
    const ypad = Math.max(20, (ymax - ymin) * .14), y0 = ymin - ypad, y1 = ymax + ypad;
    const px = (t) => left + Math.max(0, t) / Math.max(1, durationS) * (width - left - right);
    const py = (y) => top + (y1 - y) / Math.max(1, y1 - y0) * (height - top - bottom);
    const trace = (pts, color, cssClass) => {
        if (!pts.length)
            return '';
        const stride = Math.max(1, Math.ceil(pts.length / 1000));
        const subset = pts.filter((_, i) => i % stride === 0 || i === pts.length - 1);
        const path = subset.map((s, i) => `${i ? 'L' : 'M'}${px(s.tS).toFixed(1)},${py(s.h).toFixed(1)}`).join(' ');
        return `<path class="${cssClass}" d="${path}" fill="none" stroke="${color}" stroke-width="2.3" stroke-linejoin="round"/>`;
    };
    const u = project.display.altitudeUnit;
    let ticks = '';
    for (let i = 0; i <= 4; i++) {
        const y = y0 + (y1 - y0) * i / 4, yp = py(y);
        ticks += `<line x1="${left}" y1="${yp.toFixed(1)}" x2="${width - right}" y2="${yp.toFixed(1)}" stroke="#dbe4ef" stroke-width=".8"/>`;
        ticks += `<text x="${left - 8}" y="${(yp + 4).toFixed(1)}" text-anchor="end" font-size="11" fill="#607086">${Math.round(heightForDisplay(y, u))}</text>`;
    }
    // Horizontal axis is elapsed flight time, with actual time labels and no distance offset.
    const stepCount = 4;
    for (let i = 0; i <= stepCount; i++) {
        const t = durationS * i / stepCount, x = px(t);
        ticks += `<line x1="${x.toFixed(1)}" y1="${height - bottom}" x2="${x.toFixed(1)}" y2="${height - bottom + 4}" stroke="#94a3b8"/>`;
        ticks += `<text x="${x.toFixed(1)}" y="${height - 17}" font-size="10" text-anchor="middle" fill="#607086">${fmtElapsed(t)}</text>`;
    }
    const label = timingProvisional ? '未接続・制約未成立：着陸側時刻は仮配置' : '離脱からの推定経過飛行時間';
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="離脱側予想高度（Bezier接続を含む）と着陸側逆算必要高度。横軸は飛行時間。">
 <rect x="0" y="0" width="${width}" height="${height}" fill="#fff"/><text x="9" y="12" font-size="11" fill="#64748b">高度(${u}, MSL)</text>${ticks}
 ${trace(departure, '#2563eb', 'glp-departure-altitude')}${trace(arrival, '#059669', 'glp-arrival-altitude')}
 <text x="${left}" y="${height - 2}" font-size="10" fill="#64748b">${label}</text><text x="${width - right}" y="${height - 2}" font-size="10" fill="#64748b" text-anchor="end">${timingProvisional ? '計画時間 評価不可' : '計画時間 ' + fmtElapsed(durationS)}</text>
 </svg>`;
}
