export const alt = (m, u) => u === 'ft' ? m / 0.3048 : m;
export const fmtAlt = (m, u) => `${Math.round(alt(m, u)).toLocaleString('ja-JP')} ${u}`;
export function chartSvg(plan, project, width = 780, height = 200) {
    const pred = plan.forwardAll.filter(s => Number.isFinite(s.sM) && Number.isFinite(s.h));
    const req = plan.requiredAll.filter(s => Number.isFinite(s.sM) && Number.isFinite(s.h));
    const all = [...pred, ...req];
    if (!all.length)
        return '<svg viewBox="0 0 780 200"></svg>';
    const left = 48, right = 14, top = 15, bottom = 28;
    const xmin = Math.min(...all.map(v => v.sM)), xmax = Math.max(...all.map(v => v.sM)), ymin = Math.min(...all.map(v => v.h)), ymax = Math.max(...all.map(v => v.h));
    const ypad = Math.max(20, (ymax - ymin) * .14), y0 = ymin - ypad, y1 = ymax + ypad;
    const px = (x) => left + (x - xmin) / Math.max(1, xmax - xmin) * (width - left - right);
    const py = (y) => top + (y1 - y) / Math.max(1, y1 - y0) * (height - top - bottom);
    const trace = (pts, color, dash = '') => {
        if (!pts.length)
            return '';
        const stride = Math.max(1, Math.ceil(pts.length / 700));
        const subset = pts.filter((_, i) => i % stride === 0 || i === pts.length - 1);
        const path = subset.map((s, i) => `${i ? 'L' : 'M'}${px(s.sM).toFixed(1)},${py(s.h).toFixed(1)}`).join(' ');
        return `<path d="${path}" fill="none" stroke="${color}" stroke-width="2.3" ${dash ? `stroke-dasharray="${dash}"` : ''} stroke-linejoin="round"/>`;
    };
    const u = project.display.altitudeUnit;
    let ticks = '';
    for (let i = 0; i <= 4; i++) {
        const y = y0 + (y1 - y0) * i / 4, yp = py(y);
        ticks += `<line x1="${left}" y1="${yp.toFixed(1)}" x2="${width - right}" y2="${yp.toFixed(1)}" stroke="#dbe4ef" stroke-width=".8"/>`;
        ticks += `<text x="${left - 8}" y="${(yp + 4).toFixed(1)}" text-anchor="end" font-size="11" fill="#607086">${Math.round(alt(y, u))}</text>`;
    }
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="飛行高度プロファイル">
 <rect x="0" y="0" width="${width}" height="${height}" fill="#fff"/><text x="9" y="12" font-size="11" fill="#64748b">高度(${u})</text>${ticks}
 ${trace(pred, '#2563eb')}${trace(req, '#db7a26', '5 5')}
 <text x="${left}" y="${height - 6}" font-size="11" fill="#64748b">離脱</text><text x="${width - right}" y="${height - 6}" font-size="11" fill="#64748b" text-anchor="end">対地累積距離（最終） ${Math.round(xmax)} m</text>
 </svg>`;
}
export function printMapSvg(plan, width = 780, height = 420) {
    const all = plan.forwardAll.filter(s => Number.isFinite(s.e) && Number.isFinite(s.n));
    if (all.length === 0)
        return '';
    const xmin = Math.min(...all.map(p => p.e)), xmax = Math.max(...all.map(p => p.e)), ymin = Math.min(...all.map(p => p.n)), ymax = Math.max(...all.map(p => p.n));
    const margin = 40, scale = Math.min((width - margin * 2) / Math.max(50, xmax - xmin), (height - margin * 2) / Math.max(50, ymax - ymin));
    const cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2;
    const px = (e) => width / 2 + (e - cx) * scale, py = (n) => height / 2 - (n - cy) * scale;
    const group = (which, color, dashed = false) => {
        const arr = plan.forwardAll.filter(s => s.group === which);
        if (arr.length < 2)
            return '';
        const stride = Math.max(1, Math.ceil(arr.length / 600));
        const points = arr.filter((_, i) => i % stride === 0 || i === arr.length - 1);
        const path = points.map((s, i) => `${i ? 'L' : 'M'}${px(s.e).toFixed(1)},${py(s.n).toFixed(1)}`).join(' ');
        return `<path d="${path}" fill="none" stroke="${color}" stroke-width="3" ${dashed ? 'stroke-dasharray="6 5"' : ''}/>`;
    };
    const d = plan.departure.start, l = plan.landing;
    return `<svg viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="経路図（背景地図なし）">
 <rect width="${width}" height="${height}" fill="#f6f9fc"/><path d="M ${width / 2} 0 V ${height} M 0 ${height / 2} H ${width}" stroke="#cbd5e1" stroke-dasharray="3 6"/>
 ${group('departure', '#2563eb')}${group('connector', '#8b5cf6')}${group('arrival', '#059669')}${group('final', '#059669', true)}
 <circle cx="${px(d.e)}" cy="${py(d.n)}" r="5" fill="#2563eb"/><circle cx="${px(l.e)}" cy="${py(l.n)}" r="5" fill="#059669"/>
 <text x="${px(d.e) + 7}" y="${py(d.n) - 8}" font-size="12" fill="#1e293b">離脱</text><text x="${px(l.e) + 7}" y="${py(l.n) - 8}" font-size="12" fill="#1e293b">着陸</text>
 <path d="M ${width - 35} 42 V 16 M ${width - 39} 22 L ${width - 35} 16 L ${width - 31} 22" stroke="#334155" stroke-width="2" fill="none"/><text x="${width - 39}" y="12" font-size="12">N</text>
 <text x="12" y="${height - 10}" fill="#64748b" font-size="11">ENU座標による経路図（背景地図なし・縮尺は図枠に合わせて調整）</text></svg>`;
}
