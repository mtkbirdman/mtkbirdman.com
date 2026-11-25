import type { FlightLeg, TurnLeg, PlanResult, ProjectV1, Side } from '../model/types.js';
import { norm } from '../geometry/geo.js';
/** The input field shows the value currently in force without freezing inheritance. */
export function effectiveLegSpeedMps(project: ProjectV1, leg: FlightLeg): number {
    return leg.speedOverrideMps ?? project.defaultSpeedMps;
}
/** A null direction remains inherited in the JSON. Use the actual solved segment,
 * which also handles the reversed arrival-leg design order and wind corrections.
 * Direction is only used for straight legs (turns are defined by signed angle).
 */
export function effectiveStraightDirectionDeg(project: ProjectV1, plan: PlanResult, side: Side, leg: FlightLeg): number | null {
    if (leg.type !== 'straight')
        return null;
    if (leg.directionDeg !== null)
        return norm(leg.directionDeg);
    const segments = side === 'departure' ? plan.departure.segments : plan.arrival.segments;
    const segment = segments.find(s => s.leg.id === leg.id);
    const point = segment?.samples[0];
    if (point) {
        const direction = leg.guidance === 'heading' ? point.headingDeg : point.trackDeg;
        return Number.isFinite(direction) ? norm(direction) : null;
    }
    // Graceful fallback if a preceding leg failed and no simulated segment exists.
    if (side === 'departure')
        return leg.guidance === 'heading'
            ? norm(project.departure.headingDeg) : norm(plan.departure.end.trackDeg);
    return leg.guidance === 'heading'
        ? norm(plan.finalStart.headingDeg) : norm(plan.finalStart.trackDeg);
}

/** Angle actually integrated and the true bearing at the intended turn boundary.
 * The arrival boundary is before the turn in real flight (segments are stored forward). */
export function effectiveTurnValues(plan: PlanResult, side: Side, leg: TurnLeg): {angleDeg: number; boundaryDeg: number | null} {
    const path = side === 'departure' ? plan.departure : plan.arrival;
    const segment = path.segments.find(s => s.leg.id === leg.id);
    const end = side === 'departure' ? segment?.samples.at(-1) : segment?.samples[0];
    const failed = path.notices.some(n => n.legId === leg.id && n.level === 'error');
    return {
        angleDeg: segment?.leg.type === 'turn' ? segment.leg.angleDeg : leg.angleDeg,
        boundaryDeg: end && !failed ? norm(leg.guidance === 'heading' ? end.headingDeg : end.trackDeg) : null
    };
}
