import test from 'node:test';
import assert from 'node:assert/strict';
import { createDemo, makeStraight, makeTurn, trafficPreset } from '../../assets/js/model/defaults.js';
import { buildPlan, rememberResolvedTurnAngles } from '../../assets/js/planning/plan.js';
import { LocalFrame, geographicDistance, headingVector, signed } from '../../assets/js/geometry/geo.js';
import { airGround, simulatePath, simulateLeg, headingForTrack, minTurnGround, groundSpeedValid } from '../../assets/js/flight/core.js';
import { checkedSink, createSinkModel } from '../../assets/js/flight/sink.js';
import { evaluateConnector, chooseConnector } from '../../assets/js/planning/connector.js';
import { selfIntersects, curveInRange } from '../../assets/js/geometry/bezier.js';
import { loadProject, parseProject, serializeProject } from '../../assets/js/io/project.js';
import { altitudeTimeSeries, chartSvg } from '../../assets/js/ui/chart.js';
import { checkpointTable } from '../../assets/js/ui/checkpoints.js';
import { reportHtml } from '../../assets/js/ui/print.js';
const near=(a,b,tol=1e-7)=>assert.ok(Number.isFinite(a)&&Math.abs(a-b)<=tol,`${a} != ${b}, tol ${tol}`);
const errors=r=>r.notices.filter(n=>n.level==='error');
const code=(r,c)=>assert.ok(r.notices.some(n=>n.code===c),JSON.stringify(r.notices));
const wx=(fromDeg=270,speedMps=0)=>({fromDeg,speedMps});
const context=(wind=wx(),V=25)=>({speedMps:V,wind,performance:{bestLD:34,bestLDAtMps:25}});
const state=(headingDeg=0,speedMps=25,wind=wx(),e=0,n=0)=>{const g=airGround(headingDeg,speedMps,wind);return {e,n,h:1000,headingDeg,trackDeg:g.track,bankDeg:0,speedMps,groundMps:g.speed,tS:0,sM:0};};
const connector=(a=100,b=100)=>({type:'bezier_connector',handleStartM:a,handleEndM:b,curvaturePolicy:'match-neighbors',adjustmentMode:'manual'});
function aligned(p, distance=350) {
    p.connector=null;const b=buildPlan(p).arrival.start,t=headingVector(b.trackDeg),f=new LocalFrame(p.landing.point);
    p.departure={...f.toGeo({e:b.e-distance*t.e,n:b.n-distance*t.n}),altitudeMSLM:1000,headingDeg:b.headingDeg};
    p.departureLegs=[];p.connector=connector();return p;
}
function validPreset(){const p=createDemo();p.arrivalLegs=trafficPreset('left');return aligned(p);}

test('T01–04: final uses one speed and wind solution; impossible final is not timed as zero',()=>{
    for(const wind of [wx(),wx(95,30)]) {
        const p=createDemo();p.wind=wind;p.arrivalLegs=[{...makeStraight(300),speedOverrideMps:35}];
        const r=buildPlan(p),G=35-wind.speedMps;
        near(r.finalStart.speedMps,35);near(r.finalStart.groundMps,G);near(r.arrival.end.speedMps,35);
        near(r.finalTarget.at(-1).tS-r.finalTarget[0].tS,p.landing.final.groundLengthM/G);
        near(airGround(r.finalStart.headingDeg,35,wind).track,95);
    }
    // Full audit wind condition with every flown section at 35 m/s, while base is 25.
    const full=createDemo();full.wind=wx(95,30);full.arrivalLegs=[{...makeStraight(300),speedOverrideMps:35}];
    const boundary=buildPlan(full).arrival.start,t=headingVector(95),f=new LocalFrame(full.landing.point);
    full.departure={...f.toGeo({e:boundary.e-450*t.e,n:boundary.n-450*t.n}),headingDeg:95,altitudeMSLM:1000};
    full.departureLegs=[{...makeStraight(100),speedOverrideMps:35}];full.connector=connector();
    const fullResult=buildPlan(full);assert.equal(fullResult.status,'valid',JSON.stringify(fullResult.notices));
    const predicted=fullResult.forwardAll.filter(s=>s.group==='final');assert.ok(predicted.length>1);near(predicted.at(-1).tS-predicted[0].tS,120);
    const w=35**3/(2*34*25**2)+25**2/(2*34*35);near(predicted[0].h-predicted.at(-1).h,120*w);
    for(const sample of predicted){near(sample.speedMps,35);near(sample.groundMps,5);near(sample.trackDeg,95);}
    const p=createDemo();near(buildPlan(p).finalStart.speedMps,25);p.arrivalLegs=[makeStraight()];near(buildPlan(p).finalStart.speedMps,25);
    p.wind=wx(95,24.9995);let r=buildPlan(p);code(r,'GROUND_SPEED_TOO_LOW');assert.equal(r.arrival.completed,false);assert.ok(r.finalTarget.every(s=>!Number.isFinite(s.tS)));assert.equal(r.metrics.totalTimeS,null);
    p.wind=wx(95,30);r=buildPlan(p);code(r,'WIND_TRACK_INFEASIBLE');assert.equal(r.arrival.completed,false);
    p.wind=wx();p.performance={bestLD:2,bestLDAtMps:5};p.defaultSpeedMps=100;r=buildPlan(p);code(r,'SINK_MODEL_OUT_OF_RANGE');assert.equal(r.metrics.totalTimeS,null);assert.equal(r.forwardAll.filter(s=>s.group==='final').length,0);
});

test('T05–09: arrival canonical samples, same positions/times/losses and strict final boundary',()=>{
    const p=createDemo();p.arrivalLegs=[{...makeStraight(400),directionDeg:0,checkpointKind:'downwind-start'},makeTurn('right',90),{...makeStraight(300),directionDeg:95}];
    aligned(p);let r=buildPlan(p);assert.equal(r.status,'valid',JSON.stringify(r.notices));code(r,'BOUNDARY_DISCONTINUITY');
    const raw=r.arrival.samples,pred=r.forwardAll.filter(s=>s.group==='arrival'),req=r.requiredAll.filter(s=>s.group==='arrival');
    assert.equal(raw.length,pred.length);assert.equal(raw.length,req.length);
    for(let i=0;i<raw.length;i++){
        for(const k of ['e','n','headingDeg','trackDeg','speedMps','groundMps','bankDeg']){assert.equal(pred[i][k],raw[i][k]);assert.equal(req[i][k],raw[i][k]);}
        near(pred[i].tS,req[i].tS);near(pred[i].sM,req[i].sM);near(pred[i].h-req[i].h,r.metrics.connectionMarginM);
        if(i){assert.ok(pred[i].tS>=pred[i-1].tS);assert.ok(pred[i].sM>=pred[i-1].sM);near(pred[i].tS-pred[i-1].tS,raw[i].tS-raw[i-1].tS);}
    }
    for(const seg of r.arrival.segments){assert.deepEqual(seg.start,seg.samples[0]);assert.deepEqual(seg.end,seg.samples.at(-1));}
    near(r.checkpoints.find(c=>c.checkpointKind==='downwind-start').marginM,r.metrics.connectionMarginM);
    p.arrivalLegs.at(-1).directionDeg=94;r=buildPlan(p);code(r,'BOUNDARY_STATE_INVALID');assert.equal(r.connectorReady,false);assert.equal(r.metrics.connectionMarginM,null);
});

test('T06–08: canonical arrival under wind and speed steps, both turn directions',()=>{
    for(const side of ['left','right']){
        const p=createDemo();p.wind=wx(240,5);p.arrivalLegs=trafficPreset(side);p.arrivalLegs[1].speedOverrideMps=30;p.arrivalLegs[2].speedOverrideMps=28;
        aligned(p);const r=buildPlan(p);assert.equal(r.status,'valid',JSON.stringify(r.notices));code(r,'SPEED_STEP_IDEALIZED');
        const req=r.requiredAll.filter(s=>s.group==='arrival'),pred=r.forwardAll.filter(s=>s.group==='arrival');
        assert.equal(req.length,pred.length);for(let i=0;i<req.length;i++){near(pred[i].h-req[i].h,r.metrics.connectionMarginM);assert.equal(pred[i].tS,req[i].tS);assert.equal(pred[i].headingDeg,req[i].headingDeg);}
        for(let i=1;i<r.arrival.segments.length;i++){const a=r.arrival.segments[i-1].end,b=r.arrival.segments[i].start;for(const k of ['e','n','h','sM','tS'])near(a[k],b[k]);}
    }
});

test('T10–13,43: no margin or total on incomplete/invalid route, negative margin is a warning',()=>{
    let p=validPreset(),r=buildPlan(p);assert.equal(r.status,'valid',JSON.stringify(r.notices));
    p.maxBankDeg=40;p.arrivalLegs[1].bankDeg=45;aligned(p);r=buildPlan(p);assert.equal(r.connector.feasible,true);code(r,'BANK_LIMIT_EXCEEDED');assert.equal(r.status,'invalid');assert.deepEqual(r.metrics,{totalDistanceM:null,totalTimeS:null,connectionMarginM:null});assert.ok(r.checkpoints.every(c=>c.marginM===undefined));assert.match(checkpointTable(r,p),/評価不可/);assert.match(chartSvg(r,p),/計画時間 評価不可/);
    p=validPreset();p.departure.altitudeMSLM=0;r=buildPlan(p);assert.equal(r.status,'valid');assert.ok(r.metrics.connectionMarginM<0);code(r,'ALTITUDE_DEFICIT');
    p.connector=null;r=buildPlan(p);assert.equal(r.status,'invalid');assert.equal(r.forwardAll.filter(s=>s.group==='arrival').length,0);
    r=buildPlan(createDemo());assert.equal(r.status,'unplanned');assert.equal(r.metrics.totalTimeS,null);
    p=createDemo();p.wind=wx(270,30);const turn={...makeTurn('right',450),angleInputMode:'target',targetDirectionDeg:45};p.departureLegs=[turn,makeStraight()];r=buildPlan(p);rememberResolvedTurnAngles(p,r);assert.equal(turn.angleDeg,450);assert.equal(r.departure.completed,false);assert.equal(r.departure.segments.length,0);
});

test('T14–16: unique positive wind solution only, strong heading-mode drift remains available',()=>{
    for(const W of [25,30]) {const wind=wx(270,W);assert.equal(headingForTrack(90,25,wind),null);const r=simulateLeg(state(90,25,wind),makeStraight(300),context(wind),1,'departure');assert.equal(r.completed,false);code(r,'WIND_TRACK_INFEASIBLE');const h=simulateLeg(state(90,25,wind),{...makeStraight(300),guidance:'heading'},context(wind),1,'departure');assert.equal(h.completed,true);near(h.end.e,300);}
    const strong=wx(270,30);const circle=simulateLeg(state(90,25,strong),{...makeTurn('right',360),guidance:'heading'},context(strong),1,'departure');assert.equal(circle.completed,true,JSON.stringify(circle.notices));
    const wind=wx(270,30),p=createDemo();p.wind=wind;const a=state(270,25,wind),b=state(270,25,wind,350,0);assert.equal(evaluateConnector(a,b,p,100,100).feasible,false);
    p.wind=wx(270,5);const aa=state(90,25,p.wind),bb=state(90,30,p.wind,350,0);const c=evaluateConnector(aa,bb,p,100,100);assert.equal(c.feasible,true,JSON.stringify(c.notices));code(c,'SPEED_STEP_IDEALIZED');
    bb.headingDeg=270;assert.equal(evaluateConnector(aa,bb,p,100,100).feasible,false);
});

test('T17–19: analytical minimum ground speed includes interior and reverse/wrap extremes',()=>{
    const wind=wx(0,24.999),a=359.875,b=360.125;
    assert.ok(airGround(a,25,wind).speed>.05);assert.ok(airGround(b,25,wind).speed>.05);
    near(minTurnGround(25,wind,a,b),.001,1e-12);near(minTurnGround(25,wind,b,a),.001,1e-12);near(minTurnGround(25,wind,719.875,1080.125),.001,1e-12);
    const r=simulateLeg(state(a,25,wind),makeTurn('right',89),context(wind),1,'departure');assert.equal(r.completed,false);code(r,'GROUND_SPEED_TOO_LOW');
    for(const g of [.049999,.05,.050001]){assert.equal(groundSpeedValid(g),g>=.05);const w=wx(0,25-g);assert.equal(!!headingForTrack(0,25,w),g>=.05);}
});

test('T20: independent constant-wind turn position/time/loss, left/right, reverse, heading/track, multi-turn',()=>{
    let maxPosition=0,maxTime=0,maxDistance=0,maxLoss=0,maxTarget=0;
    const rad=x=>x*Math.PI/180;
    for(const wind of [wx(),wx(270,8),wx(135,24)])for(const direction of ['left','right'])for(const rev of [1,-1])for(const guidance of ['heading','ground-track'])for(const angle of [110,450]){
        const V=25,start=state(35,V,wind),leg={...makeTurn(direction,angle),guidance},r=simulateLeg(start,leg,context(wind),rev,'departure');assert.equal(r.completed,true,JSON.stringify(r.notices));
        const omega=9.80665*Math.tan(rad(30))*(direction==='right'?1:-1)/V;
        const target=guidance==='heading'?start.headingDeg:start.trackDeg;
        const sign=(direction==='right'?1:-1)*rev;
        const expectedTarget=target+sign*angle;
        const actual=guidance==='heading'?r.end.headingDeg:r.end.trackDeg;
        maxTarget=Math.max(maxTarget,Math.abs(signed(actual-expectedTarget)));
        const h0=rad(start.headingDeg),dt=r.end.tS,h1=h0+omega*dt;
        // Independently recover cumulative terminal heading from the wind triangle.
        let expectedAngle=angle;
        if(guidance==='ground-track'){
            const t=rad(expectedTarget),we=-wind.speedMps*Math.sin(rad(wind.fromDeg)),wn=-wind.speedMps*Math.cos(rad(wind.fromDeg));
            const along=we*Math.sin(t)+wn*Math.cos(t),cross=we*Math.cos(t)-wn*Math.sin(t),G=along+Math.sqrt(V*V-cross*cross);
            const h= Math.atan2(G*Math.sin(t)-we,G*Math.cos(t)-wn)*180/Math.PI;
            expectedAngle=Math.floor(angle/360)*360+((sign*(h-start.headingDeg)%360)+360)%360;
        }
        const expectedTime=rev*rad(expectedAngle)/Math.abs(omega);
        maxTime=Math.max(maxTime,Math.abs(dt-expectedTime));
        const we=-wind.speedMps*Math.sin(rad(wind.fromDeg)),wn=-wind.speedMps*Math.cos(rad(wind.fromDeg));
        const e=V/omega*(Math.cos(h0)-Math.cos(h1))+we*dt,n=V/omega*(Math.sin(h1)-Math.sin(h0))+wn*dt;
        maxPosition=Math.max(maxPosition,Math.hypot(r.end.e-e,r.end.n-n));
        const w=V/(2*34)*(1+1/Math.cos(rad(30))**2);maxLoss=Math.max(maxLoss,Math.abs(start.h-r.end.h-w*dt));
        let sum=0,N=12800;for(let i=0;i<=N;i++){const psi=h0+omega*dt*i/N;sum+=(i===0||i===N?1:i%2?4:2)*Math.hypot(V*Math.sin(psi)+we,V*Math.cos(psi)+wn);}maxDistance=Math.max(maxDistance,Math.abs(r.end.sM-sum*dt/N/3));
    }
    assert.ok(maxPosition<.01&&maxTime<.001&&maxLoss<.01&&maxDistance<.01&&maxTarget<.02,JSON.stringify({maxPosition,maxTime,maxDistance,maxLoss,maxTarget}));
    console.log('independent turn maximum errors',JSON.stringify({maxPosition,maxTime,maxDistance,maxLoss,maxTarget}));
});

// Independent de Casteljau evaluation, derivative polygons and Simpson integration.
function independentBezier(a,b,p,hA,hB,N=12800){
    const tA=headingVector(a.trackDeg),tB=headingVector(b.trackDeg),A=[a.e,a.n],B=[b.e,b.n];
    const cp=[A,[a.e+hA*tA.e,a.n+hA*tA.n],[a.e+2*hA*tA.e,a.n+2*hA*tA.n],[b.e-2*hB*tB.e,b.n-2*hB*tB.n],[b.e-hB*tB.e,b.n-hB*tB.n],B];
    const diff=(x,k)=>x.slice(1).map((v,i)=>v.map((y,j)=>k*(y-x[i][j]))),d=diff(cp,5),dd=diff(d,4);
    const dc=(points,u)=>{let row=points;while(row.length>1)row=row.slice(1).map((v,i)=>v.map((x,j)=>(1-u)*row[i][j]+u*x));return row[0];};
    const V=a.speedMps,W=p.wind.speedMps,wr=p.wind.fromDeg*Math.PI/180,we=-W*Math.sin(wr),wn=-W*Math.cos(wr);let L=0,T=0,H=0,max=0;
    for(let i=0;i<=N;i++){const u=i/N,[x,y]=dc(d,u),[xx,yy]=dc(dd,u),ds=Math.hypot(x,y),k=(x*yy-y*xx)/ds**3,along=(x*we+y*wn)/ds,cross=(-y*we+x*wn)/ds,q=Math.sqrt(V*V-cross*cross),G=along+q,bank=Math.atan(-k*G*G/q*V/9.80665),w=V**3/(2*p.performance.bestLD*p.performance.bestLDAtMps**2)+p.performance.bestLDAtMps**2/(2*p.performance.bestLD*V*Math.cos(bank)**2),coef=i===0||i===N?1:i%2?4:2;L+=coef*ds;T+=coef*ds/G;H+=coef*ds/G*w;max=Math.max(max,Math.abs(bank)*180/Math.PI);}
    return {length:L/(3*N),time:T/(3*N),loss:H/(3*N),max};
}
test('T21–23: independent Bezier convergence and between-sample bank limit',()=>{
    const p=createDemo();p.wind=wx(270,8);const a=state(15,25,p.wind),b=state(80,25,p.wind,900,1000);
    const ref=independentBezier(a,b,p,450,350,12800),fine=independentBezier(a,b,p,450,350,25600);near(ref.length,fine.length,1e-7);near(ref.loss,fine.loss,1e-7);
    const r=evaluateConnector(a,b,p,450,350);assert.equal(r.feasible,true,JSON.stringify(r.notices));near(r.lengthM,fine.length,.01);near(r.end.tS,fine.time,.001);near(r.lossM,fine.loss,.01);near(r.maxBankDeg,fine.max,.001);
    const sampleMax=Math.max(...r.samples.map(s=>Math.abs(s.bankDeg)));assert.ok(fine.max>sampleMax);p.maxBankDeg=(fine.max+sampleMax)/2;const rejected=evaluateConnector(a,b,p,450,350);assert.equal(rejected.feasible,false);code(rejected,'BANK_LIMIT_EXCEEDED');
    p.maxBankDeg=60;const degenerate=evaluateConnector(a,{...a},p,1,1);assert.equal(degenerate.feasible,false);
    console.log('independent Bezier errors',JSON.stringify({length:r.lengthM-fine.length,time:r.end.tS-fine.time,loss:r.lossM-fine.loss,maxBank:r.maxBankDeg-fine.max,limit:p.maxBankDeg}));
});

test('T24–26: open polyline crossing, touching, collinear overlap and ordinary shared endpoints',()=>{
    const pts=a=>a.map(([e,n])=>({e,n}));
    for(const x of [[[0,0],[1,1],[0,1],[1,0]],[[0,0],[1,0],[1,1],[0,0]],[[0,0],[2,0],[1,0]],[[0,0],[2,0],[2,1],[1,0],[3,0]],[[0,0],[0,0],[1,0]]])assert.equal(selfIntersects(pts(x)),true,JSON.stringify(x));
    for(const x of [[[0,0],[1,0],[2,0]],[[0,0],[1,0],[1,1],[2,1]]])assert.equal(selfIntersects(pts(x)),false);
});

test('T27–29: schema enums are strict; model-invalid inputs remain editable',()=>{
    const p=createDemo();p.departureLegs=[makeTurn()];p.connector=connector();
    for(const mutate of [q=>q.departureLegs[0].guidance=['heading'],q=>q.departureLegs[0].direction=['left'],q=>q.connector.adjustmentMode=['auto'],q=>q.display.altitudeUnit=['ft'],q=>delete q.departureLegs[0].directionDeg,q=>q.departureLegs[0].directionDeg='0']){
        const q=structuredClone(p);mutate(q);assert.throws(()=>loadProject(q));assert.throws(()=>parseProject(JSON.stringify(q)));}
    for(const v of [null,0,360]){const q=structuredClone(p);q.departureLegs[0].directionDeg=v;assert.deepEqual(parseProject(serializeProject(q)),q);}
    p.engineVersion='0.1.18-test';p.wind.speedMps=30;assert.deepEqual(loadProject(p),p);assert.equal(buildPlan(p).status,'invalid');p.schemaVersion='1.2';assert.throws(()=>loadProject(p));
});

test('T30–33: checked sink applicability independent of mathematical formula and propagates invalid altitude',()=>{
    near(createSinkModel({bestLD:2,bestLDAtMps:5}).rate(100,0),10000.0625);
    for(const [ratio,expected] of [[.099,'normal'],[.1,'normal'],[.10001,'warning'],[.2,'warning'],[.20001,'error'],[NaN,'error'],[Infinity,'error']]){
        const n=[],w=checkedSink({rate:()=>ratio*20},20,0,n);assert.equal(n[0]?.level??'normal',expected);assert.equal(w===null,expected==='error');}
    for(const bank of [0,15,30,45,60]){const n=[];assert.ok(checkedSink(createSinkModel({bestLD:34,bestLDAtMps:25}),25,bank,n)>0);assert.equal(n.length,0);}
    const p=createDemo();p.performance={bestLD:2,bestLDAtMps:5};p.defaultSpeedMps=100;
    const ctx={speedMps:100,wind:wx(),performance:p.performance},s=state(0,100);
    for(const leg of [makeStraight(300),makeTurn('right',90)]){const r=simulateLeg(s,leg,ctx,1,'departure');assert.equal(r.completed,true);code(r,'SINK_MODEL_OUT_OF_RANGE');assert.equal(r.samples[0].altitudeValid,false);assert.equal(r.end.altitudeValid,false);}
    const c=evaluateConnector(s,state(0,100,wx(),0,350),p,100,100);assert.equal(c.feasible,false);assert.equal(c.altitudeValid,false);
    const r=buildPlan(p);assert.equal(r.status,'invalid');code(r,'SINK_MODEL_OUT_OF_RANGE');
});

test('T34–38: geodesic precheck, dateline, generated-path envelope, cyclic angle',async()=>{
    const p=createDemo();p.landing.point={latDeg:0,lonDeg:0};p.departure={...p.departure,latDeg:0,lonDeg:180};p.departureLegs=[makeStraight()];let r=buildPlan(p);code(r,'LOCAL_FRAME_OUT_OF_RANGE');assert.equal(r.departureInRange,false);assert.equal(r.departure.samples.length,0);assert.ok(!r.checkpoints.some(c=>c.group==='departure'));
    const origin={latDeg:0,lonDeg:179.9};near(geographicDistance(origin,{latDeg:0,lonDeg:-179.9}),22239.016,0.001);
    for(const d of [49999,50001]){p.landing.point={latDeg:0,lonDeg:0};p.departure={...p.departure,latDeg:0,lonDeg:d/6371008.8*180/Math.PI};assert.equal(buildPlan(p).departureInRange,d<50000);}
    let s=state(90,25,wx(),49990,0);let leg=makeStraight(100);r=simulateLeg(s,leg,context(),1,'departure');code(r,'LOCAL_FRAME_OUT_OF_RANGE');assert.equal(r.completed,false);
    r=simulateLeg(s,makeTurn('right',90),context(),1,'departure');code(r,'LOCAL_FRAME_OUT_OF_RANGE');assert.equal(r.completed,false);
    assert.equal(curveInRange([{e:49990,n:0},{e:50100,n:0},{e:49990,n:0}]),false);
    const out=evaluateConnector(state(90,25,wx(),49900,0),state(270,25,wx(),49900,100),createDemo(),500,500);assert.equal(out.feasible,false);code(out,'LOCAL_FRAME_OUT_OF_RANGE');
    const f=new LocalFrame({latDeg:90,lonDeg:0});for(const xy of [{e:0,n:0},{e:100,n:100}]){const q=f.toGeo(xy);assert.ok(Number.isFinite(q.latDeg)&&Number.isFinite(q.lonDeg));}
    // Exercise map projection at the pole and preserve the original out-of-range anchor.
    const oldL=globalThis.L,paths=[],markers=[],fits=[];
    const finite=(v)=>{assert.ok(v.flat(Infinity).every(Number.isFinite));};
    const layer=(pts,opts)=>{finite(pts);return {addTo(){return this;},on(){return this;},bindTooltip(){return this;}};};
    globalThis.L={map(){return {setView(){return this;},on(){return this;},invalidateSize(){},fitBounds(v){fits.push(v);}};},svg(){return {};},tileLayer(){return {addTo(){}};},layerGroup(){return {addTo(){return this;},clearLayers(){}};},polyline(pts,opts){paths.push({pts,opts});return layer(pts,opts);},circle:layer,circleMarker:layer,marker(pts,opts){markers.push(pts);return layer(pts,opts);},divIcon(v){return v;},latLngBounds(pts){finite(pts);return pts;}};
    try {
        const {PlannerMap}=await import('../../assets/js/ui/map.js');const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){}});
        const q=createDemo();q.landing.point={latDeg:0,lonDeg:0};q.departure={...q.departure,latDeg:0,lonDeg:180};let z=buildPlan(q);map.redraw(q,z);map.fit(q,z);
        assert.ok(markers.some(v=>v[0]===0&&v[1]===180));assert.ok(!paths.some(v=>v.opts.color==='#2563eb'));assert.ok(!fits.at(-1).some(v=>v[1]===180));
        q.landing.point={latDeg:90,lonDeg:0};q.departure={...q.departure,latDeg:89.999,lonDeg:0};q.wind=wx(270,3);z=buildPlan(q);map.redraw(q,z);map.fit(q,z);
    } finally {globalThis.L=oldL;}
    near(Math.abs(signed(359.9-.1)),.2);near(signed(360),0);near(Math.abs(signed(180)),180);
    const u=headingVector(.1),wrap=evaluateConnector({...state(359.9),trackDeg:.1},state(.1,25,wx(),350*u.e,350*u.n),createDemo(),100,100);assert.equal(wrap.feasible,true,JSON.stringify(wrap.notices));
});

test('T39: coarse winner fails detailed limit; next candidate is evaluated and accepted',()=>{
    const p=createDemo();p.wind=wx(162.56983478553593,2.6410288140177727);p.maxBankDeg=10.797274307942327;
    const a=state(159.58405821584165,25,p.wind),b=state(327.16614838689566,25,p.wind,1914.5258085336536,1925.551439053379);
    const d=Math.hypot(b.e,b.n),h1=d*.16,h2=d*.27;
    const coarse=evaluateConnector(a,b,p,h1,h2,true),fine=evaluateConnector(a,b,p,h1,h2,false);
    assert.equal(coarse.feasible,true);assert.equal(fine.feasible,false);code(fine,'BANK_LIMIT_EXCEEDED');
    const candidates=[];for(const x of [.16,.27,.4,.58,.8,1.1])for(const y of [.16,.27,.4,.58,.8,1.1]){const r=evaluateConnector(a,b,p,d*x,d*y,true);candidates.push({x,y,score:(r.feasible?0:1e7)+120*r.maxBankDeg+r.lengthM+Math.max(0,r.maxBankDeg-p.maxBankDeg)*1e4});}
    candidates.sort((a,b)=>a.score-b.score);assert.equal(candidates[0].x,.16);assert.equal(candidates[0].y,.27);
    const chosen=chooseConnector(a,b,p);assert.equal(chosen.feasible,true,JSON.stringify(chosen.notices));near(chosen.handleStartM,d*.27);near(chosen.handleEndM,d*.27);
});

test('T40: deterministic bounded connector search, invalid common conditions',()=>{
    const p=createDemo();let a=state(0),b=state(0,25,wx(),0,350);const r=chooseConnector(a,b,p),s=chooseConnector(a,b,p);assert.equal(r.feasible,true);assert.equal(r.handleStartM,s.handleStartM);assert.equal(r.handleEndM,s.handleEndM);
    p.maxBankDeg=5;b=state(180,25,wx(),0,350);const z=chooseConnector(a,b,p);assert.equal(z.feasible,false);code(z,'CONNECTOR_NO_SOLUTION');
    p.wind=wx(270,30);a=state(90,25,p.wind);b=state(90,25,p.wind,350,0);assert.equal(chooseConnector(a,b,p).feasible,false);
});

test('T41–43: errors sorted before warnings; invalid chart and printed conditions are consistent',async()=>{
    const p=validPreset();p.departureLegs=Array.from({length:15},(_,i)=>({...makeStraight(1),directionDeg: i%2?90:95,speedOverrideMps:i%2?25:26}));p.maxBankDeg=5;
    const r=buildPlan(p);assert.equal(r.status,'invalid');assert.equal(r.notices[0].level,'error');let last=0;const rank={error:0,warning:1,info:2};for(const n of r.notices){assert.ok(rank[n.level]>=last);last=rank[n.level];}
    const doc=globalThis.document;globalThis.document={querySelectorAll:()=>[],querySelector:()=>null};try{
        const {App}=await import('../../assets/js/main.js');const app=Object.create(App.prototype);app.project=p;app.plan=r;app.notice='';app.updateCheckpointSelection=()=>{};const nodes=new Map();app.root={querySelector(k){if(!nodes.has(k))nodes.set(k,{innerHTML:''});return nodes.get(k);}};app.renderResults();
        const metrics=nodes.get('.glp-metrics').innerHTML,diagnostics=nodes.get('.glp-diagnostics').innerHTML;assert.equal((metrics.match(/評価不可/g)||[]).length,3);assert.match(diagnostics,/エラー \d+ 件/);assert.match(diagnostics,/件省略/);assert.doesNotMatch(metrics,/NaN|Infinity/);
        const mapEl={clientWidth:100,clientHeight:100,cloneNode(){return {querySelectorAll(){return []},removeAttribute(){},classList:{add(){}},style:{},outerHTML:'<div></div>'}}};const html=reportHtml({project:p,plan:r,mapEl,metricsHtml:metrics,diagnosticsHtml:diagnostics,title:'GLP',engineVersion:'0.1.21-test'});assert.ok(html.includes(metrics)&&html.includes(diagnostics));assert.match(html,/ファイナル実効速度/);assert.match(html,/計画時間 評価不可/);
    }finally{globalThis.document=doc;}
    const r2=buildPlan(createDemo());r2.status='invalid';r2.requiredAll=[{...r2.finalTarget[0],tS:0},{...r2.finalTarget[1],tS:1,altitudeValid:false},{...r2.finalTarget[2],tS:2}];const series=altitudeTimeSeries(r2);assert.equal(series.arrival.length,1);
});

test('T44–46: strict raw sample budget including zero, exact, short, and large forward/reverse paths',()=>{
    const ctx=context(),s=state();for(const rev of [1,-1])for(const leg of [makeStraight(20),makeTurn('right',1),{...makeTurn('left',1),guidance:'heading'}]){
        const full=simulateLeg(s,leg,ctx,rev,'departure'),n=full.samples.length;assert.equal(full.completed,true);for(const budget of [0,1,n-1,n]){const r=simulateLeg(s,leg,ctx,rev,'departure',budget);assert.ok(r.samples.length<=budget);assert.equal(r.completed,budget>=n);}}
    const r=simulatePath(s,[makeStraight(20),makeStraight(20),makeStraight(20)],ctx,'departure',false,3);assert.equal(r.completed,false);assert.equal(r.segments.length,1);code(r,'PATH_COMPLEXITY_LIMIT');
    for(const reverse of [false,true]){
        const start=performance.now(),legs=[12500,12500,4999.25,1].map(a=>({...makeTurn('right',a),guidance:'heading'}));if(reverse)legs.reverse();
        const z=simulatePath(s,legs,ctx,reverse?'arrival':'departure',reverse);const raw=z.segments.reduce((n,l)=>n+l.samples.length,0);assert.equal(z.completed,false);assert.equal(raw,120000);assert.ok(raw<=120000);code(z,'PATH_COMPLEXITY_LIMIT');console.log('budget load',JSON.stringify({reverse,raw,ms:performance.now()-start,rss:process.memoryUsage().rss}));
    }
    const tooLong=simulateLeg(s,{...makeTurn('right',100001),guidance:'heading'},ctx,1,'departure');assert.equal(tooLong.completed,false);code(tooLong,'TURN_COMPLEXITY_LIMIT');
});

test('GLP-015: regularity certificate, limits and bounded subdivision', async()=>{
    const {curveRegular,controlPoints,TANGENT_MIN_M}=await import('../../assets/js/geometry/bezier.js');
    const cp=b=>controlPoints({e:0,n:0},{e:0,n:1000},{e:0,n:1},{e:0,n:1},0,0,200,b);
    assert.equal(curveRegular(cp(580),0),'uncertain');
    assert.equal(curveRegular(cp(580),16,1),'uncertain');
    assert.equal(curveRegular(cp(580)),'regular');
    for(const b of [590.625,590.6251])assert.notEqual(curveRegular(cp(b)),'regular');
    const zero=cp(200);zero[1]={...zero[0]};assert.equal(curveRegular(zero),'below-limit');
    zero[1]={e:TANGENT_MIN_M/10,n:0};assert.equal(curveRegular(zero),'below-limit');
    zero[1].e=Infinity;assert.equal(curveRegular(zero),'uncertain');
    assert.equal(curveRegular([]),'uncertain');
    assert.equal(curveRegular(cp(200),16,0),'uncertain');
    const curved=controlPoints({e:0,n:0},{e:900,n:1000},headingVector(15),headingVector(80),.0001,-.0001,450,350);
    assert.equal(curveRegular(curved),'regular');
});

test('GLP-015: stationary tangent and tiny reversal rejected in coarse and detailed evaluation',()=>{
    const p=createDemo();
    for(const h of [0,90,237])for(const b of [590.625,590.6251])for(const coarse of [true,false]){
        const t=headingVector(h),A=state(h,25,wx(),200,300),B=state(h,25,wx(),200+1000*t.e,300+1000*t.n);
        const r=evaluateConnector(A,B,p,200,b,coarse);
        assert.equal(r.feasible,false);assert.equal(r.altitudeValid,false);assert.equal(r.samples.length,0);
        assert.ok(r.points.length>2);assert.ok(r.notices.some(n=>['CONNECTOR_TANGENT_TOO_SMALL','CONNECTOR_REGULARITY_UNCERTAIN'].includes(n.code)));
        assert.equal(r.lengthM,0);assert.equal(r.lossM,0);assert.deepEqual(r.end,A);
    }
    const a=state(0),b=state(0,25,wx(),0,1000);
    assert.equal(evaluateConnector(a,b,p,200,580).feasible,true);
    assert.equal(chooseConnector(a,b,p).feasible,true);
});

function stationaryPlan(){
    const p=createDemo();p.landing.point={latDeg:0,lonDeg:0};p.landing.runwayTrackDeg=90;
    p.arrivalLegs=[makeStraight(300)];
    p.departure={...p.departure,latDeg:0,lonDeg:Math.asin(-1900/6378137)*180/Math.PI,headingDeg:90};
    p.connector=connector(200,590.6251);return p;
}
test('GLP-015: buildPlan retains geometry only, and recovers on editing',()=>{
    const p=stationaryPlan(),r=buildPlan(p);
    assert.equal(r.connectorReady,true);assert.equal(r.status,'invalid');assert.equal(r.connector.feasible,false);
    assert.deepEqual(r.metrics,{totalDistanceM:null,totalTimeS:null,connectionMarginM:null});
    assert.ok(r.connector.points.length>2);assert.equal(r.connector.samples.length,0);
    assert.ok(r.forwardAll.every(s=>s.group!=='connector'));assert.ok(r.requiredAll.length>0);
    assert.ok(r.checkpoints.every(c=>c.marginM===undefined));
    p.connector.handleEndM=200;assert.equal(buildPlan(p).status,'valid');
});

test('GLP-016: completed checkpoints only, preserving partial segments and custom names',()=>{
    for(const reverse of [false,true]){
        const p=createDemo();p.landing.point={latDeg:0,lonDeg:0};p.landing.runwayTrackDeg=180;
        const unfinished={...makeStraight(reverse?50000:2000),checkpointLabel:'未到達の名称',...(reverse?{checkpointKind:'downwind-start'}:{})};
        const complete=makeStraight(reverse?100:100);
        if(reverse)p.arrivalLegs=[unfinished,complete];
        else {p.departure={...p.departure,latDeg:.442,lonDeg:0,headingDeg:0};p.departureLegs=[complete,unfinished];}
        const r=buildPlan(p),side=reverse?r.arrival:r.departure;
        assert.equal(r.status,'invalid');assert.equal(side.completed,false);
        const seg=side.segments.find(s=>s.leg.id===unfinished.id);assert.ok(seg);assert.equal(seg.completed,false);assert.ok(seg.samples.length>1);
        assert.ok(r.checkpoints.some(c=>c.legId===complete.id));assert.ok(!r.checkpoints.some(c=>c.legId===unfinished.id));
        assert.doesNotMatch(checkpointTable(r,p),/未到達の名称/);assert.equal(unfinished.checkpointLabel,'未到達の名称');
        assert.ok(r.checkpoints.some(c=>c.kind==='final-start'));assert.equal(r.metrics.totalTimeS,null);
    }
    const p=createDemo();p.wind=wx(0,30);p.departureLegs=[makeStraight()];
    const r=buildPlan(p);assert.equal(r.departure.completed,false);assert.ok(!r.checkpoints.some(c=>c.legId===p.departureLegs[0].id));
});

test('GLP-017: generated and fallback handles remain within the persistence contract',()=>{
    const p=createDemo();p.connector=connector();
    for(const [a,b] of [[state(0,25,wx(),0,-44999),state(0,25,wx(),0,45000)], [state(0,25,wx(),0,-49999),state(0,25,wx(),0,49999)], [state(0),state(0,25,wx(),0,.2)], [state(0),state(0,25,wx(),0,Infinity)], [state(0,25,wx(),0,-1000000),state(0,25,wx(),0,1000000)]]){
        const r=chooseConnector(a,b,p);
        for(const h of [r.handleStartM,r.handleEndM])assert.ok(Number.isFinite(h)&&h>=1&&h<=50000);
        p.connector.handleStartM=r.handleStartM;p.connector.handleEndM=r.handleEndM;
        assert.deepEqual(parseProject(serializeProject(p)),p);
        if(!Number.isFinite(b.n)||Math.abs(b.n)>50000)assert.equal(r.feasible,false);
    }
    for(const h of [0,50000.1,NaN,Infinity]){
        const r=evaluateConnector(state(0),state(0,25,wx(),0,1000),p,h,100);
        assert.equal(r.feasible,false);assert.equal(r.samples.length,0);code(r,'CONNECTOR_HANDLE_OUT_OF_RANGE');
    }
});
