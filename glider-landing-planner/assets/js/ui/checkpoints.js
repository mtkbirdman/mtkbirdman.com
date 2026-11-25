import { fmtAlt } from './chart.js?v=0.1.16';
import { escapeHtml } from './escape-html.js?v=0.1.16';
export function landingRelativeHeight(checkpoint, project) {
    return checkpoint.altitudeM - project.landing.elevationMSLM;
}
export function checkpointRows(plan, project) {
    const rows = [];
    for (const cp of plan.checkpoints) {
        if (cp.checkpointKind === 'downwind-start')
            rows.push({ kind: 'margin', label: 'ダウンウィンド開始の高度余裕', relativeM: cp.marginM, legId: cp.legId });
        rows.push({ kind: 'height', label: cp.label, relativeM: cp.altitudeValid === false ? undefined : landingRelativeHeight(cp, project), mslM: cp.altitudeValid === false ? undefined : cp.altitudeM, legId: cp.legId });
    }
    return rows;
}
export function checkpointTable(plan, project, print = false) {
    const unit = project.display.altitudeUnit;
    const tableClass = print ? 'glp-print-table' : 'glp-table';
    const body = checkpointRows(plan, project).map(r => {
        const relative = r.relativeM === undefined ? '評価不可' : r.kind === 'margin' && r.relativeM > 0 ? `+${fmtAlt(r.relativeM, unit)}` : fmtAlt(r.relativeM, unit);
        const msl = r.mslM === undefined ? '—' : fmtAlt(r.mslM, unit);
        const classes = [r.kind === 'margin' ? 'glp-cp-margin-row' : ''].filter(Boolean).join(' ');
        const style = `${classes ? ` class="${classes}"` : ''}${r.legId ? ` data-leg-id="${escapeHtml(r.legId)}"` : ''}`;
        return `<tr${style}><td>${escapeHtml(r.label)}</td><td>${relative}</td><td>${msl}</td></tr>`;
    }).join('');
    return `<table class="${tableClass}"><thead><tr><th>地点</th><th>滑空場基準高度（${unit}）</th><th>MSL高度（${unit}）</th></tr></thead><tbody>${body}</tbody></table>`;
}
