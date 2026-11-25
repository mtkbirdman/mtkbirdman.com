import { heightForDisplay } from './height-units.js?v=0.1.14';
export const fmtAlt = (m, u) => `${Math.round(heightForDisplay(m, u)).toLocaleString('ja-JP')} ${u}`;
export const fmtElapsed = (seconds) => {
    const s = Math.max(0, Math.round(seconds));
    return `${Math.floor(s / 60)}:${String(s % 60).padStart(2, '0')}`;
};
/**
 * Use genuine elapsed seconds, not ENU distance or a distance-based display offset.
 * Reverse-integrated arrival samples carry negative tS; shift them to the actual
 * end of the connector, without modifying the canonical flight/plan data.
 * The two curves intentionally remain separate even when their x coordinates meet.
 */
export function altitudeTimeSeries(plan) {
    const finite = (s) => s.altitudeValid !== false && Number.isFinite(s.tS) && Number.isFinite(s.h);
    const departure = [...plan.departure.samples, ...(plan.connector?.samples.slice(1) ?? [])].filter(finite);
    const joinedTime = departure.at(-1)?.tS ?? 0;
    const arrivalStartRaw = plan.arrival.samples[0]?.tS ?? plan.arrival.start.tS;
    const arrivalShift = joinedTime - arrivalStartRaw;
    const arrival = plan.arrival.samples.filter(finite).map(s => ({ ...s, tS: s.tS + arrivalShift }));
    const arrivalEnd = arrival.at(-1)?.tS ?? joinedTime;
    const finalStartRaw = plan.finalTarget[0]?.tS ?? 0;
    arrival.push(...plan.finalTarget.slice(1).filter(finite).map(s => ({ ...s, tS: arrivalEnd + s.tS - finalStartRaw })));
    let durationS = 0;
    for (const s of departure)
        durationS = Math.max(durationS, s.tS);
    for (const s of arrival)
        durationS = Math.max(durationS, s.tS);
    // Time and geometry can remain evaluable even when the polar gives no altitude.
    for (const s of plan.forwardAll)
        if (Number.isFinite(s.tS))
            durationS = Math.max(durationS, s.tS);
    return { departure, arrival, durationS, timingProvisional: !plan.connector || !plan.connector.feasible };
}
export function chartSvg(plan, project, width = 780, height = 200) {
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
 <text x="${left}" y="${height - 2}" font-size="10" fill="#64748b">${label}</text><text x="${width - right}" y="${height - 2}" font-size="10" fill="#64748b" text-anchor="end">計画時間 ${fmtElapsed(durationS)}</text>
 </svg>`;
}
