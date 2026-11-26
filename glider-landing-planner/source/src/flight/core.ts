import type { WindConfig, PerformanceConfig, State, FlightLeg, TurnLeg, Sample, Notice, SegmentResult, PathResult, Point2 } from '../model/types.js';
import { createSinkModel, checkedSink } from './sink.js';
import type { SinkModel } from './sink.js';
import { rad, deg, norm, signed, headingVector, courseOf, magnitude, dot, mul, add, inLocalRange, LOCAL_RADIUS_M } from '../geometry/geo.js';
export const MIN_GROUND_MPS = 0.05;
// Tiny tolerance only for cancellation at the inclusive decimal boundary.
export const groundSpeedValid = (g:number):boolean => Number.isFinite(g) && g >= MIN_GROUND_MPS-1e-12;
export const PATH_SAMPLE_BUDGET = 120000;
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
    if (!(wind.speedMps < speedMps)) return null;
    const t = headingVector(trackDeg), w = windVector(wind), parallel = dot(w, t), perp2 = dot(w, w) - parallel ** 2;
    const q2 = speedMps ** 2 - Math.max(0, perp2);
    if (!(q2 > 0))
        return null;
    const groundMps = parallel + Math.sqrt(q2);
    if (!groundSpeedValid(groundMps))
        return null;
    const va = { e: groundMps * t.e - w.e, n: groundMps * t.n - w.n };
    return { headingDeg: courseOf(va.e, va.n), groundMps };
}
/** Convenience wrapper for independent numerical checks; plan integration reuses one evaluator. */
export function sinkRate(speedMps: number, bankDeg: number, p: PerformanceConfig): number {
    return createSinkModel(p).rate(speedMps, bankDeg);
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
/** Exact minimum G over a monotone unwrapped heading interval. */
export function minTurnGround(V: number, wind: WindConfig, a: number, b: number): number {
    const lo = Math.min(a,b), hi = Math.max(a,b);
    const extreme = wind.fromDeg + 360*Math.ceil((lo-wind.fromDeg)/360);
    return extreme <= hi ? Math.abs(V-wind.speedMps) : Math.min(airGround(a,V,wind).speed,airGround(b,V,wind).speed);
}
function turnInRange(s: State, dt: number, bank: number, ctx: FlightContext, extra:{remaining:number}, depth=0): boolean {
    if (!inLocalRange(s)) return false;
    if (magnitude(s)+(s.speedMps+ctx.wind.speedMps)*Math.abs(dt) <= LOCAL_RADIUS_M) return true;
    if(extra.remaining--<=0)return false;
    const end=advanceTurn(s,dt,bank,ctx,null);
    if (!inLocalRange(end)) return false;
    if (Math.min(magnitude(s),magnitude(end))+(s.speedMps+ctx.wind.speedMps)*Math.abs(dt) <= LOCAL_RADIUS_M) return true;
    if (depth>=12) return false;
    const mid=advanceTurn(s,dt/2,bank,ctx,null);
    return turnInRange(s,dt/2,bank,ctx,extra,depth+1)&&turnInRange(mid,dt/2,bank,ctx,extra,depth+1);
}
export function simulateLeg(initial: State, leg: FlightLeg, ctx: FlightContext, timeDirection: 1 | -1, group: Sample['group'], budget=PATH_SAMPLE_BUDGET): {
    samples: Sample[]; end: State; notices: Notice[]; effectiveLeg: FlightLeg; completed: boolean;
} {
    const notices: Notice[] = [], samples: Sample[] = [];
    const V=leg.speedOverrideMps??ctx.speedMps;
    let current=reState({...initial,speedMps:V},ctx.wind);
    const effectiveLeg:FlightLeg=leg.type==='turn'&&leg.angleInputMode==='target' ? {...leg,angleDeg:closestTurnAngle(leg.targetDirectionDeg!,leg.guidance==='heading'?current.headingDeg:current.trackDeg,leg.direction,timeDirection===-1,leg.angleDeg)} : leg;
    const result=(completed:boolean)=>({samples,end:current,notices,effectiveLeg,completed});
    const fail=(code:string,message:string)=>{notices.push({code,level:'error',message,legId:leg.id});return result(false);};
    const budgetFail=()=>fail('PATH_COMPLEXITY_LIMIT','経路の計算点数が上限に達しました。');
    const rangeFail=()=>fail('LOCAL_FRAME_OUT_OF_RANGE','経路が局所計算領域（50 km）を超えるか、領域内であることを確認できません。');
    const groundFail=()=>fail('GROUND_SPEED_TOO_LOW','対地速度が0.05 m/s未満です。');
    if (budget<1) return budgetFail();
    if (!inLocalRange(current)) return rangeFail();
    if (leg.guidance==='ground-track' && !(ctx.wind.speedMps<V)) return fail('WIND_TRACK_INFEASIBLE','航跡方位指定では風速が対気速度未満である必要があります。');
    if (Math.abs(initial.speedMps-V)>1e-6) notices.push({code:'SPEED_STEP_IDEALIZED',level:'warning',message:'速度変更は瞬間的と仮定（加減速・エネルギー変化を未評価）。',legId:leg.id});
    const model=ctx.sinkModel??createSinkModel(ctx.performance);
    const validG=groundSpeedValid;
    if (leg.type==='straight') {
        const target=leg.directionDeg??(leg.guidance==='heading'?current.headingDeg:current.trackDeg);
        const desired=leg.guidance==='heading'?target:headingForTrack(target,V,ctx.wind)?.headingDeg;
        if (desired===undefined) return groundFail();
        if (Math.abs(signed(desired-current.headingDeg))>0.3) notices.push({code:'BOUNDARY_DISCONTINUITY',level:'warning',message:'直線境界の機首方位変化は瞬間的と仮定します。',legId:leg.id});
        current=reState({...current,headingDeg:norm(desired),bankDeg:0},ctx.wind);
        const g=airGround(current.headingDeg,V,ctx.wind);
        if (!validG(g.speed)) return groundFail();
        const steps=Math.max(1,Math.min(2000,Math.ceil(leg.lengthM/20)));
        if (steps+1>budget) return budgetFail();
        const dt=timeDirection*leg.lengthM/g.speed/steps;
        const sink=checkedSink(model,V,0,notices,leg.id);
        if (sink===null) current={...current,altitudeValid:false};
        samples.push(makeSample(current,leg.id,group));
        for(let i=0;i<steps;i++) {
            const next={...current,e:current.e+g.e*dt,n:current.n+g.n*dt,h:current.h-(sink!==null&&current.altitudeValid!==false?sink*dt:0),tS:current.tS+dt,sM:current.sM+g.speed*dt};
            if (!inLocalRange(next)) return rangeFail();
            current=next;samples.push(makeSample(current,leg.id,group));
        }
        return result(true);
    }
    if (!validG(current.groundMps)) return groundFail();
    const turn=effectiveLeg as TurnLeg, bank=leg.bankDeg*(leg.direction==='right'?1:-1), sign=(leg.direction==='right'?1:-1)*timeDirection;
    const omega=GRAVITY*Math.tan(rad(bank))/V;
    if (!Number.isFinite(omega)||Math.abs(omega)<1e-7) return fail('TURN_INVALID','旋回率が不正です。');
    let angle=turn.angleDeg;
    if (leg.guidance==='ground-track') {
        const target=norm(current.trackDeg+sign*angle), sol=headingForTrack(target,V,ctx.wind);
        if (!sol) return groundFail();
        const turns=Math.floor(angle/360), remainder=angle%360;
        angle=turns*360+(remainder<1e-10?0:norm(sign*(sol.headingDeg-current.headingDeg)));
    }
    if (!Number.isFinite(angle)||angle<=0||angle>TURN_STEP_BUDGET*TURN_MAX_STEP_DEG) return fail('TURN_COMPLEXITY_LIMIT','旋回角が計算可能なステップ数を超えました。');
    const step=Math.max(.25,angle/TURN_STEP_BUDGET), steps=Math.min(TURN_STEP_BUDGET,Math.max(1,Math.ceil(angle/step)));
    if (steps+1>budget) return budgetFail();
    if (!validG(minTurnGround(V,ctx.wind,current.headingDeg,current.headingDeg+sign*angle))) return groundFail();
    const sink=checkedSink(model,V,bank,notices,leg.id);
    current={...current,bankDeg:bank,...(sink===null?{altitudeValid:false}:{})};
    samples.push(makeSample(current,leg.id,group));
    const dt=timeDirection*rad(angle)/Math.abs(omega)/steps;
    let progress=0;
    const rangeBudget={remaining:20000};
    for(let i=0;i<steps;i++) {
        if (!turnInRange(current,dt,bank,ctx,rangeBudget)) return rangeFail();
        const next=advanceTurn(current,dt,bank,ctx,sink);
        if (!inLocalRange(next)) return rangeFail();
        // W<V makes the track monotone, but a near-singular step can exceed 180°.
        const delta=sign*(leg.guidance==='heading'?next.headingDeg-current.headingDeg:next.trackDeg-current.trackDeg);
        progress+=Math.abs(delta)<1e-10?0:norm(delta);
        current=next;samples.push(makeSample(current,leg.id,group));
    }
    if (Math.abs(progress-turn.angleDeg)>.02) return fail('TURN_TARGET_NOT_REACHED','指定した旋回角に到達しませんでした。');
    return result(true);
}
export function simulatePath(start: State, legs: FlightLeg[], ctx: FlightContext, group: 'departure'|'arrival', reverse=false, budget=PATH_SAMPLE_BUDGET): PathResult {
    let current={...start}, remaining=Math.min(PATH_SAMPLE_BUDGET,Math.max(0,budget)), completed=true;
    ctx={...ctx,sinkModel:ctx.sinkModel??createSinkModel(ctx.performance)};
    const notices:Notice[]=[], segments:SegmentResult[]=[];
    for(const leg of reverse?[...legs].reverse():legs) {
        const r=simulateLeg(current,leg,ctx,reverse?-1:1,group,remaining);
        current=r.end;remaining-=r.samples.length;notices.push(...r.notices);
        const samples=reverse?[...r.samples].reverse():r.samples;
        if(samples.length) {
            const seg={leg:r.effectiveLeg,samples,start:{...samples[0]},end:{...samples.at(-1)!},completed:r.completed};
            if(reverse) segments.unshift(seg);else segments.push(seg);
        }
        if(!r.completed){completed=false;break;}
    }
    for(let i=1;i<segments.length;i++) {
        const a=segments[i-1].end,b=segments[i].start;
        if(Math.hypot(a.e-b.e,a.n-b.n)>1e-6 || Math.abs(a.h-b.h)>1e-6 || Math.abs(a.tS-b.tS)>1e-6 || Math.abs(a.sM-b.sM)>1e-6) {
            completed=false;notices.push({code:'BOUNDARY_STATE_INVALID',level:'error',message:'隣接レグの位置・高度・時刻・距離が不連続です。',legId:segments[i].leg.id});
        }
        if(Math.abs(signed(a.headingDeg-b.headingDeg))>.3 && !notices.some(n=>n.code==='BOUNDARY_DISCONTINUITY'&&n.legId===segments[i].leg.id)) notices.push({code:'BOUNDARY_DISCONTINUITY',level:'warning',message:'隣接レグの機首方位変化は瞬間的と仮定します。',legId:segments[i].leg.id});
    }
    const samples=segments.flatMap((s,i)=>i?s.samples.slice(1):s.samples);
    return {samples:samples.length?samples:legs.length?[]:[makeSample(start,group+'-start',group)],segments,start:segments[0]?.start??start,end:segments.at(-1)?.end??start,notices,completed};
}
