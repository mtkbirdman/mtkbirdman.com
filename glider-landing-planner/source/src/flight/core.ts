import type { WindConfig, PerformanceConfig, State, FlightLeg, TurnLeg, Sample, Notice, SegmentResult, PathResult, Point2 } from '../model/types.js';
import { createSinkModel } from './sink.js';
import type { SinkModel } from './sink.js';
import { rad, deg, norm, signed, headingVector, courseOf, magnitude, dot, mul, add } from '../geometry/geo.js';
export const GRAVITY = 9.80665;
/** No fixed angle limit in the project format; bound computation instead. */
export const TURN_STEP_BUDGET = 50000;
export const TURN_MAX_STEP_DEG = 2;
export interface FlightContext {
    speedMps: number;
    wind: WindConfig;
    performance: PerformanceConfig;
    sinkModel?: SinkModel;
}
export function windVector(w: WindConfig): Point2 { return mul(headingVector(w.fromDeg), -w.speedMps); }
export function airGround(headingDeg: number, speedMps: number, wind: WindConfig): {
    e: number;
    n: number;
    speed: number;
    track: number;
} {
    const w = windVector(wind), v = add(mul(headingVector(headingDeg), speedMps), w);
    return { ...v, speed: magnitude(v), track: courseOf(v.e, v.n) };
}
export function headingForTrack(trackDeg: number, speedMps: number, wind: WindConfig): {
    headingDeg: number;
    groundMps: number;
} | null {
    const t = headingVector(trackDeg), w = windVector(wind), parallel = dot(w, t), perp2 = dot(w, w) - parallel ** 2;
    const q2 = speedMps ** 2 - Math.max(0, perp2);
    if (q2 <= 1e-8)
        return null;
    const groundMps = parallel + Math.sqrt(q2);
    if (groundMps <= 1e-7)
        return null;
    const va = { e: groundMps * t.e - w.e, n: groundMps * t.n - w.n };
    return { headingDeg: courseOf(va.e, va.n), groundMps };
}
/** Convenience wrapper for independent numerical checks; plan integration reuses one evaluator. */
export function sinkRate(speedMps: number, bankDeg: number, p: PerformanceConfig): number {
    return createSinkModel(p).rate(speedMps, bankDeg);
}
function safeSink(model: SinkModel, speedMps: number, bankDeg: number, notices: Notice[], legId?: string): number | null {
    const rate = model.rate(speedMps, bankDeg);
    if (Number.isFinite(rate) && rate > 0) return rate;
    notices.push({code: 'SINK_RESULT_INVALID', level: 'error', message: '沈下率が正の有限値になりません。', ...(legId ? {legId} : {})});
    return null;
}
/** Choose the nearest positive number of turns that ends at the desired true bearing.
 * Arrival legs are designed backwards, but direction still means real-flight left/right.
 * In an exact tie, retain the smaller positive solution. */
export function closestTurnAngle(targetDeg: number, boundaryDeg: number, direction: 'left' | 'right', reverse: boolean, previousAngleDeg: number): number {
    const sign = direction === 'right' ? 1 : -1;
    const remainder = norm(sign * (reverse ? boundaryDeg - targetDeg : targetDeg - boundaryDeg));
    const base = remainder < 1e-9 || 360 - remainder < 1e-9 ? 360 : remainder;
    const k = Math.max(0, Math.floor((previousAngleDeg - base) / 360));
    const lower = base + 360 * k;
    const upper = lower + 360;
    return Math.abs(previousAngleDeg - lower) <= Math.abs(previousAngleDeg - upper) ? lower : upper;
}
export function reState(s: State, wind: WindConfig): State {
    const g = airGround(s.headingDeg, s.speedMps, wind);
    return { ...s, trackDeg: g.track, groundMps: g.speed };
}
export function makeSample(s: State, segmentId: string, group: Sample['group']): Sample { return { ...s, segmentId, group }; }
function advanceTurn(s: State, dt: number, bankSigned: number, ctx: FlightContext, sink: number | null): State {
    const V = s.speedMps, w = windVector(ctx.wind), omega = GRAVITY * Math.tan(rad(bankSigned)) / V;
    // RK4 position integration for constant bank and wind, stable in both time directions.
    const velocity = (heading: number) => ({ e: V * Math.sin(rad(heading)) + w.e, n: V * Math.cos(rad(heading)) + w.n });
    const k1 = velocity(s.headingDeg), k2 = velocity(s.headingDeg + deg(omega * dt / 2)), k4 = velocity(s.headingDeg + deg(omega * dt));
    const e = s.e + dt * (k1.e + 4 * k2.e + k4.e) / 6;
    const n = s.n + dt * (k1.n + 4 * k2.n + k4.n) / 6;
    const h = sink === null || s.altitudeValid === false ? s.h : s.h - sink * dt;
    const out = reState({ ...s, e, n, h, headingDeg: norm(s.headingDeg + deg(omega * dt)), tS: s.tS + dt, bankDeg: bankSigned, ...(sink === null || s.altitudeValid === false ? { altitudeValid: false } : {}) }, ctx.wind);
    const avgG = (s.groundMps + out.groundMps) / 2;
    out.sM = s.sM + avgG * dt;
    return out;
}
export function simulateLeg(initial: State, leg: FlightLeg, ctx: FlightContext, timeDirection: 1 | -1, group: Sample['group']): {
    samples: Sample[];
    end: State;
    notices: Notice[];
    effectiveLeg: FlightLeg;
} {
    const notices: Notice[] = [];
    const V = leg.speedOverrideMps ?? ctx.speedMps;
    let current = reState({ ...initial, speedMps: V }, ctx.wind);
    // Resolve the saved design intent at the actual leg speed, before turn integration.
    // This matters for ground-track guidance in a crosswind.
    const effectiveLeg: FlightLeg = leg.type === 'turn' && leg.angleInputMode === 'target'
        ? { ...leg, angleDeg: closestTurnAngle(leg.targetDirectionDeg!,
            leg.guidance === 'heading' ? current.headingDeg : current.trackDeg,
            leg.direction, timeDirection === -1, leg.angleDeg) }
        : leg;
    const sinkModel = ctx.sinkModel ?? createSinkModel(ctx.performance);
    if (Math.abs(initial.speedMps - V) > 1e-6) {
        notices.push({ code: 'SPEED_STEP_IDEALIZED', level: 'warning', message: '速度変更は瞬間的と仮定（加減速・エネルギー変化を未評価）。', legId: leg.id });
    }
    if (leg.type === 'straight') {
        const target = leg.directionDeg ?? (leg.guidance === 'heading' ? current.headingDeg : current.trackDeg);
        const desiredHeading = leg.guidance === 'heading' ? target : headingForTrack(target, V, ctx.wind)?.headingDeg;
        if (desiredHeading === undefined) {
            notices.push({ code: 'WIND_TRACK_INFEASIBLE', level: 'error', message: 'この風・速度では指定した対地直線航跡を維持できません。', legId: leg.id });
            return { samples: [makeSample(current, leg.id, group)], end: current, notices, effectiveLeg };
        }
        if (Math.abs(signed(desiredHeading - current.headingDeg)) > 0.1) {
            notices.push({ code: 'BOUNDARY_DISCONTINUITY', level: 'warning', message: `直線開始境界で機首方位が ${Math.abs(signed(desiredHeading - current.headingDeg)).toFixed(1)}° 変化します（遷移未モデル化）。`, legId: leg.id });
        }
        current = reState({ ...current, headingDeg: norm(desiredHeading), bankDeg: 0 }, ctx.wind);
        const samples = [makeSample(current, leg.id, group)];
        const g = airGround(current.headingDeg, V, ctx.wind);
        if (g.speed < 1e-6) {
            notices.push({ code: 'GROUND_SPEED_ZERO', level: 'error', message: '対地速度がゼロのため直線距離を積算できません。', legId: leg.id });
            return { samples, end: current, notices, effectiveLeg };
        }
        const sink = safeSink(sinkModel, V, 0, notices, leg.id);
        if (sink === null) current = { ...current, altitudeValid: false };
        const steps = Math.max(1, Math.min(2000, Math.ceil(leg.lengthM / 20)));
        const dt = timeDirection * leg.lengthM / g.speed / steps;
        for (let i = 0; i < steps; i++) {
            current = { ...current, e: current.e + g.e * dt, n: current.n + g.n * dt, h: current.h - (sink !== null && current.altitudeValid !== false ? sink * dt : 0), tS: current.tS + dt, sM: current.sM + g.speed * dt };
            samples.push(makeSample(current, leg.id, group));
        }
        return { samples, end: current, notices, effectiveLeg };
    }
    const turnLeg = effectiveLeg as TurnLeg;
    // Long thermalling turns can exceed any number of revolutions. Reject only runs whose
    // required integration resolution would exceed the browser CPU/memory budget.
    const stepDeg = Math.max(0.25, turnLeg.angleDeg / TURN_STEP_BUDGET);
    if (stepDeg > TURN_MAX_STEP_DEG) {
        notices.push({ code: 'TURN_COMPLEXITY_LIMIT', level: 'error', message: '旋回角が計算可能なステップ数を超えました。レグを分割してください。', legId: leg.id });
        return { samples: [makeSample(current, leg.id, group)], end: current, notices, effectiveLeg };
    }
    const bankSigned = leg.bankDeg * (leg.direction === 'right' ? 1 : -1);
    const omega = GRAVITY * Math.tan(rad(bankSigned)) / V;
    if (!Number.isFinite(omega) || Math.abs(omega) < 1e-7) {
        notices.push({ code: 'TURN_INVALID', level: 'error', message: 'バンク角が小さすぎて旋回を計算できません。', legId: leg.id });
        return { samples: [makeSample(current, leg.id, group)], end: current, notices, effectiveLeg };
    }
    const turnSink = safeSink(sinkModel, V, Math.abs(bankSigned), notices, leg.id);
    current = { ...current, bankDeg: bankSigned, ...(turnSink === null ? { altitudeValid: false } : {}) };
    const samples = [makeSample(current, leg.id, group)];
    if (leg.guidance === 'heading') {
        const steps = Math.max(1, Math.ceil(turnLeg.angleDeg / stepDeg));
        const dt = timeDirection * rad(turnLeg.angleDeg) / Math.abs(omega) / steps;
        for (let i = 0; i < steps; i++) {
            current = advanceTurn(current, dt, bankSigned, ctx, turnSink);
            samples.push(makeSample(current, leg.id, group));
        }
    }
    else {
        const signWanted = (leg.direction === 'right' ? 1 : -1) * timeDirection;
        let turned = 0;
        const stepDt = timeDirection * rad(stepDeg) / Math.abs(omega);
        for (let i = 0; i < TURN_STEP_BUDGET + 1; i++) {
            const next = advanceTurn(current, stepDt, bankSigned, ctx, turnSink);
            if (next.groundMps < 0.05) {
                notices.push({ code: 'TURN_TARGET_NOT_REACHED', level: 'error', message: '対地速度がほぼゼロになりました。', legId: leg.id });
                break;
            }
            const delta = signed(next.trackDeg - current.trackDeg);
            const progress = signWanted * delta;
            if (progress < -0.05) {
                notices.push({ code: 'TURN_TARGET_NOT_REACHED', level: 'error', message: '風により対地航跡角の変化が非単調です。', legId: leg.id });
                break;
            }
            if (turned + progress >= turnLeg.angleDeg) {
                let lo = 0, hi = 1;
                for (let j = 0; j < 25; j++) {
                    const mid = (lo + hi) / 2, q = advanceTurn(current, stepDt * mid, bankSigned, ctx, turnSink);
                    const p = signWanted * signed(q.trackDeg - current.trackDeg);
                    if (turned + p < turnLeg.angleDeg)
                        lo = mid;
                    else
                        hi = mid;
                }
                current = advanceTurn(current, stepDt * hi, bankSigned, ctx, turnSink);
                samples.push(makeSample(current, leg.id, group));
                turned = turnLeg.angleDeg;
                break;
            }
            turned += progress;
            current = next;
            samples.push(makeSample(current, leg.id, group));
        }
        if (turned < turnLeg.angleDeg - 0.02) {
            notices.push({ code: 'TURN_TARGET_NOT_REACHED', level: 'error', message: '指定した対地航跡角の変化量に到達しませんでした（風条件または計算ステップ数の制約）。', legId: leg.id });
        }
    }
    return { samples, end: current, notices, effectiveLeg };
}
export function simulatePath(start: State, legs: FlightLeg[], ctx: FlightContext, group: 'departure' | 'arrival', reverse = false): PathResult {
    let current = { ...start };
    ctx = { ...ctx, sinkModel: ctx.sinkModel ?? createSinkModel(ctx.performance) };
    const notices: Notice[] = [];
    const segments: SegmentResult[] = [];
    const order = reverse ? [...legs].reverse() : legs;
    let sampleBudget = 0;
    for (const leg of order) {
        if (sampleBudget >= 120000) {
            notices.push({ code: 'PATH_COMPLEXITY_LIMIT', level: 'error', message: '経路の計算点数が上限に達しました。レグ数や旋回角を減らしてください。' });
            break;
        }
        const begin = { ...current };
        const r = simulateLeg(current, leg, ctx, reverse ? -1 : 1, group);
        current = r.end;
        sampleBudget += r.samples.length;
        notices.push(...r.notices);
        if (reverse)
            segments.unshift({ leg: r.effectiveLeg, samples: [...r.samples].reverse(), start: { ...current }, end: begin });
        else
            segments.push({ leg: r.effectiveLeg, samples: r.samples, start: begin, end: { ...current } });
    }
    const samples = segments.flatMap((segment, index) => index === 0 ? segment.samples : segment.samples.slice(1));
    return { samples: samples.length ? samples : [makeSample(start, group + '-start', group)], segments, start: reverse ? current : start, end: reverse ? start : current, notices };
}
