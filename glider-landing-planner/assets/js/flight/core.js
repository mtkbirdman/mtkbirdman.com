import { rad, deg, norm, signed, headingVector, courseOf, magnitude, dot, mul, add } from '../geometry/geo.js';
export const GRAVITY = 9.80665;
export function windVector(w) { return mul(headingVector(w.fromDeg), -w.speedMps); }
export function airGround(headingDeg, speedMps, wind) {
    const w = windVector(wind), v = add(mul(headingVector(headingDeg), speedMps), w);
    return { ...v, speed: magnitude(v), track: courseOf(v.e, v.n) };
}
export function headingForTrack(trackDeg, speedMps, wind) {
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
export function sinkRate(speedMps, bankDeg, p) {
    if (p.model === 'fixed-ld')
        return speedMps / p.fixedLD;
    const Vbest = p.bestLDAtMps, E = p.bestLD;
    const a = 1 / (2 * E * Vbest ** 2), b = Vbest ** 2 / (2 * E);
    return a * speedMps ** 3 + b / (speedMps * Math.cos(rad(bankDeg)) ** 2);
}
export function reState(s, wind) {
    const g = airGround(s.headingDeg, s.speedMps, wind);
    return { ...s, trackDeg: g.track, groundMps: g.speed };
}
export function makeSample(s, segmentId, group) { return { ...s, segmentId, group }; }
function advanceTurn(s, dt, bankSigned, ctx) {
    const V = s.speedMps, w = windVector(ctx.wind), omega = GRAVITY * Math.tan(rad(bankSigned)) / V;
    // RK4 position integration for constant bank and wind, stable in both time directions.
    const velocity = (heading) => ({ e: V * Math.sin(rad(heading)) + w.e, n: V * Math.cos(rad(heading)) + w.n });
    const k1 = velocity(s.headingDeg), k2 = velocity(s.headingDeg + deg(omega * dt / 2)), k4 = velocity(s.headingDeg + deg(omega * dt));
    const e = s.e + dt * (k1.e + 4 * k2.e + k4.e) / 6;
    const n = s.n + dt * (k1.n + 4 * k2.n + k4.n) / 6;
    const h = s.h - sinkRate(V, Math.abs(bankSigned), ctx.performance) * dt;
    const out = reState({ ...s, e, n, h, headingDeg: norm(s.headingDeg + deg(omega * dt)), tS: s.tS + dt, bankDeg: bankSigned }, ctx.wind);
    const avgG = (s.groundMps + out.groundMps) / 2;
    out.sM = s.sM + avgG * dt;
    return out;
}
export function simulateLeg(initial, leg, ctx, timeDirection, group) {
    const notices = [];
    const V = leg.speedOverrideMps ?? ctx.speedMps;
    let current = reState({ ...initial, speedMps: V }, ctx.wind);
    if (Math.abs(initial.speedMps - V) > 1e-6) {
        notices.push({ code: 'SPEED_STEP_IDEALIZED', level: 'warning', message: '速度変更は瞬間的と仮定（加減速・エネルギー変化を未評価）。', legId: leg.id });
    }
    if (leg.type === 'straight') {
        const target = leg.directionDeg ?? (leg.guidance === 'heading' ? current.headingDeg : current.trackDeg);
        const desiredHeading = leg.guidance === 'heading' ? target : headingForTrack(target, V, ctx.wind)?.headingDeg;
        if (desiredHeading === undefined) {
            notices.push({ code: 'WIND_TRACK_INFEASIBLE', level: 'error', message: 'この風・速度では指定した対地直線航跡を維持できません。', legId: leg.id });
            return { samples: [makeSample(current, leg.id, group)], end: current, notices };
        }
        if (Math.abs(signed(desiredHeading - current.headingDeg)) > 0.1) {
            notices.push({ code: 'BOUNDARY_DISCONTINUITY', level: 'warning', message: `直線開始境界で機首方位が ${Math.abs(signed(desiredHeading - current.headingDeg)).toFixed(1)}° 変化します（遷移未モデル化）。`, legId: leg.id });
        }
        current = reState({ ...current, headingDeg: norm(desiredHeading), bankDeg: 0 }, ctx.wind);
        const samples = [makeSample(current, leg.id, group)];
        const g = airGround(current.headingDeg, V, ctx.wind);
        if (g.speed < 1e-6) {
            notices.push({ code: 'GROUND_SPEED_ZERO', level: 'error', message: '対地速度がゼロのため直線距離を積算できません。', legId: leg.id });
            return { samples, end: current, notices };
        }
        const steps = Math.max(1, Math.min(2000, Math.ceil(leg.lengthM / 20)));
        const dt = timeDirection * leg.lengthM / g.speed / steps;
        for (let i = 0; i < steps; i++) {
            current = { ...current, e: current.e + g.e * dt, n: current.n + g.n * dt, h: current.h - sinkRate(V, 0, ctx.performance) * dt, tS: current.tS + dt, sM: current.sM + g.speed * dt };
            samples.push(makeSample(current, leg.id, group));
        }
        return { samples, end: current, notices };
    }
    const bankSigned = leg.bankDeg * (leg.direction === 'right' ? 1 : -1);
    const omega = GRAVITY * Math.tan(rad(bankSigned)) / V;
    if (!Number.isFinite(omega) || Math.abs(omega) < 1e-7) {
        notices.push({ code: 'TURN_INVALID', level: 'error', message: 'バンク角が小さすぎて旋回を計算できません。', legId: leg.id });
        return { samples: [makeSample(current, leg.id, group)], end: current, notices };
    }
    current = { ...current, bankDeg: bankSigned };
    const samples = [makeSample(current, leg.id, group)];
    if (leg.guidance === 'heading') {
        const steps = Math.max(1, Math.ceil(leg.angleDeg / 0.25));
        const dt = timeDirection * rad(leg.angleDeg) / Math.abs(omega) / steps;
        for (let i = 0; i < steps; i++) {
            current = advanceTurn(current, dt, bankSigned, ctx);
            samples.push(makeSample(current, leg.id, group));
        }
    }
    else {
        const signWanted = (leg.direction === 'right' ? 1 : -1) * timeDirection;
        let turned = 0;
        const stepDt = timeDirection * rad(0.25) / Math.abs(omega);
        for (let i = 0; i < Math.ceil(1080 / 0.25); i++) {
            const next = advanceTurn(current, stepDt, bankSigned, ctx);
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
            if (turned + progress >= leg.angleDeg) {
                let lo = 0, hi = 1;
                for (let j = 0; j < 25; j++) {
                    const mid = (lo + hi) / 2, q = advanceTurn(current, stepDt * mid, bankSigned, ctx);
                    const p = signWanted * signed(q.trackDeg - current.trackDeg);
                    if (turned + p < leg.angleDeg)
                        lo = mid;
                    else
                        hi = mid;
                }
                current = advanceTurn(current, stepDt * hi, bankSigned, ctx);
                samples.push(makeSample(current, leg.id, group));
                turned = leg.angleDeg;
                break;
            }
            turned += progress;
            current = next;
            samples.push(makeSample(current, leg.id, group));
        }
        if (turned < leg.angleDeg - 0.02) {
            notices.push({ code: 'TURN_TARGET_NOT_REACHED', level: 'error', message: '指定した対地航跡角の変化量に到達しませんでした。', legId: leg.id });
        }
    }
    if (ctx.performance.model === 'fixed-ld')
        notices.push({ code: 'TURN_SINK_UNMODELED', level: 'info', message: '固定L/Dモデルでは旋回による追加沈下率を評価していません。', legId: leg.id });
    return { samples, end: current, notices };
}
export function simulatePath(start, legs, ctx, group, reverse = false) {
    let current = { ...start };
    const notices = [];
    const segments = [];
    const order = reverse ? [...legs].reverse() : legs;
    for (const leg of order) {
        const begin = { ...current };
        const r = simulateLeg(current, leg, ctx, reverse ? -1 : 1, group);
        current = r.end;
        notices.push(...r.notices);
        if (reverse)
            segments.unshift({ leg, samples: [...r.samples].reverse(), start: { ...current }, end: begin });
        else
            segments.push({ leg, samples: r.samples, start: begin, end: { ...current } });
    }
    const samples = reverse ? segments.flatMap((seg, i) => i === 0 ? seg.samples : seg.samples.slice(1)) : segments.flatMap((seg, i) => i === 0 ? seg.samples : seg.samples.slice(1));
    return { samples: samples.length ? samples : [makeSample(start, group + '-start', group)], segments, start: reverse ? current : start, end: reverse ? start : current, notices };
}
