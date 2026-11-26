import { validHandle, clampHandle } from '../model/limits.js';
import type { State, Point2, ConnectionResult, Notice, Sample, ProjectV1 } from '../model/types.js';
import { headingVector, courseOf, rad, deg, dot, magnitude, mul, add, sub, leftNormal, signed, inLocalRange } from '../geometry/geo.js';
import { point, derivative, curvature, controlPoints, selfIntersects, curveInRange, curveRegular, TANGENT_MIN_M } from '../geometry/bezier.js';
import { windVector, GRAVITY, MIN_GROUND_MPS, groundSpeedValid, headingForTrack, airGround } from '../flight/core.js';
import { createSinkModel, checkedSink } from '../flight/sink.js';
const EPS = 1e-8;
const GROUND_RANGE_MAX_DEPTH = 16, GROUND_RANGE_NODE_BUDGET = 20000;
type GroundRange = 'valid' | 'below-limit' | 'uncertain';
/** Sufficient certificate on the original derivative, not a sampled minimum.
 * Called only after the curve and its endpoint wind solutions are validated. */
export function connectorGroundRange(cp: Point2[], V: number, project: ProjectV1,
    maxDepth = GROUND_RANGE_MAX_DEPTH, budget = GROUND_RANGE_NODE_BUDGET): GroundRange {
    const W = project.wind.speedMps;
    if (!(Number.isFinite(V) && Number.isFinite(W) && W >= 0 && W < V) ||
        cp.length !== 6 || !cp.every(p => Number.isFinite(p.e) && Number.isFinite(p.n))) return 'uncertain';
    // Reverse triangle inequality: all headings satisfy the existing inclusive threshold.
    if (groundSpeedValid(V - W)) return 'valid';
    if (W === 0) return 'below-limit';
    const g0 = MIN_GROUND_MPS - 1e-12;
    const unit = mul(windVector(project.wind), 1 / W);
    const c = (g0 * g0 - (V - W) * (V + W)) / (2 * g0 * W);
    if (!Number.isFinite(c) || !Number.isFinite(unit.e) || !Number.isFinite(unit.n)) return 'uncertain';
    let remaining = Math.max(0, Math.min(GROUND_RANGE_NODE_BUDGET, Math.floor(budget)));
    const depthLimit = Math.max(0, Math.min(GROUND_RANGE_MAX_DEPTH, Math.floor(maxDepth)));
    const visit = (q: Point2[], depth: number): GroundRange => {
        if (!(remaining > 0)) return 'uncertain';
        remaining--;
        const norms = q.map(magnitude), projections = q.map(v => dot(unit, v));
        if (![...norms, ...projections].every(Number.isFinite)) return 'uncertain';
        // Endpoints are actual derivative values; a failed wind solution is observed evidence.
        for (const i of [0, q.length - 1]) {
            if (!(norms[i] > TANGENT_MIN_M)) return 'uncertain';
            if (!headingForTrack(courseOf(q[i].e, q[i].n), V, project.wind)) return 'below-limit';
        }
        const axis = q.reduce((a, b) => add(a, b), {e: 0, n: 0}), len = magnitude(axis);
        const lower = Number.isFinite(len) && len > 0 ? Math.max(0, Math.min(...q.map(v => dot(v, mul(axis, 1 / len))))) : 0;
        const upper = Math.max(...norms);
        // Require W_hat·D >= c |D|. The sign of c determines which norm bound is safe.
        // Strict comparison here is conservative at floating-point boundaries.
        if (Math.min(...projections) > c * (c < 0 ? lower : upper)) return 'valid';
        if (!(depth < depthLimit)) return 'uncertain';
        const left = [q[0]], right = [q.at(-1)!]; let row = q;
        while (row.length > 1) {
            row = row.slice(1).map((v, i) => mul(add(v, row[i]), .5));
            left.push(row[0]); right.unshift(row.at(-1)!);
        }
        const a = visit(left, depth + 1); if (a === 'below-limit') return a;
        const b = visit(right, depth + 1);
        return b === 'below-limit' ? b : a === 'valid' && b === 'valid' ? 'valid' : 'uncertain';
    };
    return visit(cp.slice(1).map((v, i) => mul(sub(v, cp[i]), 5)), 0);
}
/** Conditions common to every handle candidate, independent of curve geometry. */
export function connectorBoundaryValid(start: State, end: State, project: ProjectV1): boolean {
    if (!inLocalRange(start)||!inLocalRange(end)||!groundSpeedValid(start.groundMps)||!groundSpeedValid(end.groundMps)) return false;
    const a=headingForTrack(start.trackDeg,start.speedMps,project.wind);
    const c=headingForTrack(end.trackDeg,start.speedMps,project.wind);
    const actualA=airGround(start.headingDeg,start.speedMps,project.wind);
    const actualB=airGround(end.headingDeg,end.speedMps,project.wind);
    const same=Math.abs(start.speedMps-end.speedMps)<=1e-4;
    return !!a&&!!c&&groundSpeedValid(actualA.speed)&&groundSpeedValid(actualB.speed)
        &&Math.abs(actualA.speed-start.groundMps)<=1e-4&&Math.abs(actualB.speed-end.groundMps)<=1e-4
        &&Math.abs(signed(a.headingDeg-start.headingDeg))<=.3&&Math.abs(signed(actualB.track-end.trackDeg))<=.3
        &&Math.abs(signed(c.headingDeg-end.headingDeg))<=(same ? .3 : 179.7);

}
export function endpointCurvature(state: State, bankDeg: number, project: ProjectV1): number {
    if (Math.abs(bankDeg) < 0.01)
        return 0;
    const v = state.speedMps, omega = GRAVITY * Math.tan(rad(bankDeg)) / v;
    const Vg = add(mul(headingVector(state.headingDeg), v), windVector(project.wind));
    const dV = { e: v * Math.cos(rad(state.headingDeg)) * omega, n: -v * Math.sin(rad(state.headingDeg)) * omega };
    return (Vg.e * dV.n - Vg.n * dV.e) / (magnitude(Vg) ** 3);
}
export function evaluateConnector(start: State, endBoundary: State, project: ProjectV1, a: number, b: number, skipHeavy = false, sinkModel = createSinkModel(project.performance)): ConnectionResult {
    if(!validHandle(a)||!validHandle(b)) return {points:[],samples:[],end:{...start},lengthM:0,lossM:0,maxBankDeg:0,feasible:false,altitudeValid:false,notices:[{code:'CONNECTOR_HANDLE_OUT_OF_RANGE',level:'error',message:'Bezierハンドル長が保存範囲外です。'}],handleStartM:a,handleEndM:b};
    const A: Point2 = { e: start.e, n: start.n }, B: Point2 = { e: endBoundary.e, n: endBoundary.n };
    const tA = headingVector(start.trackDeg), tB = headingVector(endBoundary.trackDeg);
    const kA = endpointCurvature(start, start.bankDeg, project), kB = endpointCurvature(endBoundary, endBoundary.bankDeg, project);
    const cp = controlPoints(A, B, tA, tB, kA, kB, a, b);
    const approximateSpan = magnitude(sub(B, A)) + a + b;
    const steps = skipHeavy ? 56 : Math.min(2000, Math.max(160, Math.ceil(approximateSpan / 3)));
    const points: Point2[] = [], samples: Sample[] = [], notices: Notice[] = [];
    let invalid = false, maxBank = 0, length = 0, loss = 0, time = 0;
    const w = windVector(project.wind), V = start.speedMps;
    let prev: {p:Point2;g:number;sink:number}|null=null;
    let altitudeValid=start.altitudeValid!==false;
    const commonValid=connectorBoundaryValid(start,endBoundary,project);
    if(!commonValid || !curveInRange(cp)) {
        const code=commonValid||!inLocalRange(start)||!inLocalRange(endBoundary)?'LOCAL_FRAME_OUT_OF_RANGE':'BOUNDARY_STATE_INVALID';
        return {points:[],samples:[],end:{...start},lengthM:0,lossM:0,maxBankDeg:0,feasible:false,altitudeValid:false,notices:[{code,level:'error',message:'接続の風条件・境界または局所領域が成立しません。'}],handleStartM:a,handleEndM:b};
    }
    const regularity=curveRegular(cp);
    if(regularity!=='regular') {
        const preview=Array.from({length:steps+1},(_,i)=>point(cp,i/steps));
        const code=regularity==='below-limit'?'CONNECTOR_TANGENT_TOO_SMALL':'CONNECTOR_REGULARITY_UNCERTAIN';
        const message=regularity==='below-limit'?'Bezier導関数の大きさが下限以下です。':'評価予算内でBezier曲線の正則性を確認できません。';
        return {points:preview.every(inLocalRange)?preview:[],samples:[],end:{...start},lengthM:0,lossM:0,maxBankDeg:0,feasible:false,altitudeValid:false,notices:[{code,level:'error',message}],handleStartM:a,handleEndM:b};
    }
    const groundRange = connectorGroundRange(cp, V, project);
    if (groundRange !== 'valid') {
        const preview = Array.from({length:steps+1}, (_,i) => point(cp,i/steps));
        const code = groundRange === 'below-limit' ? 'GROUND_SPEED_TOO_LOW' : 'CONNECTOR_GROUND_RANGE_UNCERTAIN';
        const message = groundRange === 'below-limit' ? '接続曲線の対地速度下限を満たしません。' : '評価予算内で接続曲線全体の対地速度下限を確認できません。';
        return {points:preview.every(inLocalRange)?preview:[],samples:[],end:{...start},lengthM:0,lossM:0,maxBankDeg:0,feasible:false,altitudeValid:false,notices:[{code,level:'error',message}],handleStartM:a,handleEndM:b};
    }
    let horizontalInvalid=false;
    const evaluate=(u:number)=>{
        const p=point(cp,u),d=derivative(cp,u),dsdu=magnitude(d),t=dsdu>EPS?mul(d,1/dsdu):tA,k=curvature(cp,u);
        const track=courseOf(t.e,t.n), sol=headingForTrack(track,V,project.wind);
        const Wn=dot(w,leftNormal(t)), q2=V*V-Wn*Wn;
        const G=sol?.groundMps??NaN;
        const bank=deg(Math.atan(-k*G*G/Math.sqrt(q2)*V/GRAVITY));
        if(!sol||!Number.isFinite(bank)||dsdu<=TANGENT_MIN_M||!inLocalRange(p)) {invalid=true;horizontalInvalid=true;if(!sol&&!notices.some(n=>n.code==='GROUND_SPEED_TOO_LOW'))notices.push({code:'GROUND_SPEED_TOO_LOW',level:'error',message:'接続曲線の対地速度下限を満たしません。'});}
        else maxBank=Math.max(maxBank,Math.abs(bank));
        const sink=checkedSink(sinkModel,V,bank,notices);
        if(sink===null){altitudeValid=false;invalid=true;}
        return {p,dsdu,G,bank,track,heading:sol?.headingDeg??NaN,sink:sink??0};
    };
    const nodes:ReturnType<typeof evaluate>[]=[];
    for(let i=0;i<=steps;i++) {
        const r=evaluate(i/steps); nodes.push(r);points.push(r.p);
        if(prev){const ds=magnitude(sub(r.p,prev.p));length+=ds;if(Number.isFinite(r.G)&&Number.isFinite(prev.g)&&groundSpeedValid(r.G)&&groundSpeedValid(prev.g)){const dt=ds*.5*(1/r.G+1/prev.g);time+=dt;if(altitudeValid)loss+=dt*.5*(r.sink+prev.sink);}else invalid=true;}
        samples.push({...r.p,h:start.h-loss,headingDeg:r.heading,trackDeg:r.track,bankDeg:r.bank,speedMps:V,groundMps:r.G,tS:start.tS+time,sM:start.sM+length,segmentId:'connector',group:'connector',altitudeValid});
        prev={p:r.p,g:r.G,sink:r.sink};
    }
    const extraPoints: {u:number;p:Point2}[]=[];
    // Refine local extrema and suspicious intervals, not the integration mesh.
    if(!skipHeavy&&!invalid&&maxBank<=project.maxBankDeg+1e-6){
        let budget=20000,uncertain=false;
        const extra=(u:number)=>{if(budget<=0){uncertain=true;return null;}budget--;const r=evaluate(u);extraPoints.push({u,p:r.p});return r;};
        const peak=(l:number,h:number)=>{for(let k=0;k<24&&h-l>1e-8&&!invalid;k++){const x=l+(h-l)*.38196601125,y=h-(h-l)*.38196601125,rx=extra(x),ry=extra(y);if(!rx||!ry)return;if(Math.abs(rx.bank)<Math.abs(ry.bank))l=x;else h=y;}};
        for(let i=1;i<steps&&!invalid;i++){const a0=Math.abs(nodes[i-1].bank),a1=Math.abs(nodes[i].bank),a2=Math.abs(nodes[i+1].bank);if(a1>=a0&&a1>=a2&&(a1>a0+1e-10||a1>a2+1e-10))peak((i-1)/steps,(i+1)/steps);}
        const refine=(l:number,h:number,a:ReturnType<typeof evaluate>,b:ReturnType<typeof evaluate>,depth:number)=>{
            if(invalid||uncertain||maxBank>project.maxBankDeg+1e-6)return;
            const prior=maxBank,m=extra((l+h)/2);if(!m)return;
            const bend=Math.abs(Math.abs(m.bank)-.5*(Math.abs(a.bank)+Math.abs(b.bank)));
            const near=[a,m,b].some(r=>r.G<.1||r.dsdu<1e-3||Math.abs(project.maxBankDeg-Math.abs(r.bank))<=.5||r.sink/V>=.19);
            if(bend>.05||(near&&(depth===0||maxBank-prior>1e-4||bend>1e-4))){if(depth>=8){uncertain=true;return;}refine(l,(l+h)/2,a,m,depth+1);refine((l+h)/2,h,m,b,depth+1);}
        };
        for(let i=0;i<steps&&!invalid&&!uncertain;i++)refine(i/steps,(i+1)/steps,nodes[i],nodes[i+1],0);
        if(uncertain){invalid=true;notices.push({code:'CONNECTOR_EVALUATION_UNCERTAIN',level:'error',message:'追加評価予算内で接続制約を確定できません。'});}
    }
    if(horizontalInvalid){samples.length=0;altitudeValid=false;}
    if(!altitudeValid) for(const s of samples)s.altitudeValid=false;
    // Golden-search evaluations can coincide to roundoff; do not create artificial zero-length edges.
    // Retain the original mesh check, and check the refined mesh including the extra extrema points.
    const intersectionPoints:Point2[]=[];
    if(extraPoints.length) for(const {p} of [...points.map((p,i)=>({u:i/steps,p})),...extraPoints].sort((a,b)=>a.u-b.u)) {
        if(!intersectionPoints.length||magnitude(sub(p,intersectionPoints.at(-1)!))>1e-5)intersectionPoints.push(p);
    }
    if (selfIntersects(points) || (intersectionPoints.length>1 && selfIntersects(intersectionPoints))) {
        invalid = true;
        notices.push({ code: 'CONNECTOR_SELF_INTERSECTION', level: 'error', message: '接続曲線が自己交差しています。' });
    }
    if (maxBank > project.maxBankDeg + 1e-6) {
        invalid = true;
        notices.push({ code: 'BANK_LIMIT_EXCEEDED', level: 'error', message: `Bezierの必要バンク角（最大${maxBank.toFixed(1)}°）が設定上限を超過。` });
    }
    if (invalid && !notices.some(n => n.level === 'error'))
        notices.push({ code: 'CONNECTOR_NO_SOLUTION', level: 'error', message: '曲線追従に必要な速度・曲率条件を満たしません。' });
    const end = samples.length ? { ...samples[samples.length - 1] } : { ...start };
    if (samples.length) {
        const startHeading = samples[0].headingDeg, endHeading = samples[samples.length - 1].headingDeg;
        const diffA=Math.abs(signed(startHeading-start.headingDeg)),diffB=Math.abs(signed(endHeading-endBoundary.headingDeg));
        const same=Math.abs(endBoundary.speedMps-V)<=1e-4;
        const actualEnd=airGround(endBoundary.headingDeg,endBoundary.speedMps,project.wind);
        if(diffA>.3||Math.abs(signed(samples[0].trackDeg-start.trackDeg))>.3||Math.abs(signed(end.trackDeg-endBoundary.trackDeg))>.3||(same&&diffB>.3)||(!same&&(!groundSpeedValid(actualEnd.speed)||Math.abs(signed(actualEnd.track-endBoundary.trackDeg))>.3||diffB>=179.7))){invalid=true;notices.push({code:'BOUNDARY_STATE_INVALID',level:'error',message:'Bezier境界の機首方位または航跡方位が不連続です。'});}
        else if(!same)notices.push({code:'SPEED_STEP_IDEALIZED',level:'warning',message:'Bezier終点で速度・機首方位の瞬間変更を仮定します。'});
    }
    return { points, samples, end, lengthM: length, lossM: loss, maxBankDeg: maxBank, feasible: !invalid, altitudeValid, notices, handleStartM: a, handleEndM: b };
}
export function chooseConnector(start: State, end: State, project: ProjectV1, sinkModel = createSinkModel(project.performance)): ConnectionResult {
    const length = Math.hypot(end.e - start.e, end.n - start.n);
    const fallback=(useSpan=true)=>{
        const h=useSpan&&Number.isFinite(length)?clampHandle(Math.max(10,length*.35)):null;
        const a=h??(validHandle(project.connector?.handleStartM)?project.connector!.handleStartM:70);
        const b=h??(validHandle(project.connector?.handleEndM)?project.connector!.handleEndM:70);
        const r=evaluateConnector(start,end,project,a,b,false,sinkModel);
        r.feasible=false;r.altitudeValid=false;
        if(!r.notices.some(n=>n.level==='error'))r.notices.push({code:'CONNECTOR_NO_SOLUTION',level:'error',message:'接続候補を確定できません。'});
        return r;
    };
    if(!Number.isFinite(length))return fallback();
    if (length < 1) {
        return evaluateConnector(start, end, project, 5, 5, false, sinkModel);
    }
    if(!connectorBoundaryValid(start,end,project))return fallback();
    const candidates:{a:number;b:number;score:number}[]=[];
    for(const x of [.16,.27,.40,.58,.8,1.1])for(const y of [.16,.27,.40,.58,.8,1.1]){
        const a=Math.max(10,length*x),b=Math.max(10,length*y);
        if(!validHandle(a)||!validHandle(b))continue;
        const r=evaluateConnector(start,end,project,a,b,true,sinkModel);
        const score=r.samples.length?(r.feasible?0:1e7)+r.maxBankDeg*120+r.lengthM+Math.max(0,r.maxBankDeg-project.maxBankDeg)*1e4:Infinity;
        candidates.push({a,b,score:Number.isFinite(score)?score:Infinity});
    }
    candidates.sort((a,b)=>a.score===b.score?0:a.score-b.score);
    let first:ConnectionResult|undefined;
    for(const c of candidates){const r=evaluateConnector(start,end,project,c.a,c.b,false,sinkModel);if(!first||(!first.points.length&&r.points.length))first=r;if(r.feasible)return r;}
    if(!first)return fallback(false);
    first.notices.push({code:'CONNECTOR_NO_SOLUTION',level:'error',message:'探索範囲内で制約を満たすBezier接続は見つかりませんでした。'});
    return first;
}
