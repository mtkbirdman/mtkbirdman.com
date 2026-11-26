import { norm } from '../geometry/geo.js?v=1.0.0';
/** The input field shows the value currently in force without freezing inheritance. */
export function effectiveLegSpeedMps(project, leg) {
    return leg.speedOverrideMps ?? project.defaultSpeedMps;
}
/** A null direction remains inherited in the JSON. Use the actual solved segment,
 * which also handles the reversed arrival-leg design order and wind corrections.
 * Direction is only used for straight legs (turns are defined by signed angle).
 */
export function effectiveStraightDirectionDeg(_project, plan, side, leg) {
    if (leg.type !== 'straight')
        return null;
    if (leg.directionDeg !== null)
        return norm(leg.directionDeg);
    const segments = side === 'departure' ? plan.departure.segments : plan.arrival.segments;
    const segment = segments.find(s => s.leg.id === leg.id);
    const point = segment?.completed ? segment.samples[0] : undefined;
    if (point) {
        const direction = leg.guidance === 'heading' ? point.headingDeg : point.trackDeg;
        return Number.isFinite(direction) ? norm(direction) : null;
    }
    return null;
}
/** Angle actually integrated and the true bearing at the intended turn boundary.
 * The arrival boundary is before the turn in real flight (segments are stored forward). */
export function effectiveTurnValues(plan, side, leg) {
    const path = side === 'departure' ? plan.departure : plan.arrival;
    const segment = path.segments.find(s => s.leg.id === leg.id);
    const end = side === 'departure' ? segment?.samples.at(-1) : segment?.samples[0];
    const failed = path.notices.some(n => n.legId === leg.id && n.level === 'error');
    const completed = !!segment?.completed && !failed;
    return {
        completed,
        angleDeg: completed && segment?.leg.type === 'turn' ? segment.leg.angleDeg : leg.angleDeg,
        boundaryDeg: end && segment?.completed && !failed ? norm(leg.guidance === 'heading' ? end.headingDeg : end.trackDeg) : null
    };
}
