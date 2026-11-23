import { fmtAlt } from './chart.js';
const encode = (value) => value.replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c] ?? c));
export function landingRelativeHeight(checkpoint, project) {
    return checkpoint.altitudeM - project.landing.elevationMSLM;
}
export function checkpointRows(plan, project) {
    const rows = [];
    for (const cp of plan.checkpoints) {
        if (cp.checkpointKind === 'downwind-start')
            rows.push({ kind: 'margin', label: 'ダウンウィンド開始の高度余裕', relativeM: cp.marginM });
        rows.push({ kind: 'height', label: cp.label, relativeM: landingRelativeHeight(cp, project), mslM: cp.altitudeM });
    }
    return rows;
}
export function checkpointTable(plan, project, print = false) {
    const unit = project.display.altitudeUnit;
    const tableClass = print ? 'glp-print-table' : 'glp-table';
    const body = checkpointRows(plan, project).map(r => {
        const relative = r.relativeM === undefined ? '評価不可' : r.kind === 'margin' && r.relativeM > 0 ? `+${fmtAlt(r.relativeM, unit)}` : fmtAlt(r.relativeM, unit);
        const msl = r.mslM === undefined ? '—' : fmtAlt(r.mslM, unit);
        const style = r.kind === 'margin' ? ' class="glp-cp-margin-row"' : '';
        return `<tr${style}><td>${encode(r.label)}</td><td>${relative}</td><td>${msl}</td></tr>`;
    }).join('');
    return `<table class="${tableClass}"><thead><tr><th>地点</th><th>滑空場基準高度（${unit}）</th><th>MSL高度（${unit}）</th></tr></thead><tbody>${body}</tbody></table>`;
}
