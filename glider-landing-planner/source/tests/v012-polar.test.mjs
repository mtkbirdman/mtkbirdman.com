import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo, makeTurn, trafficPreset} from '../../assets/js/model/defaults.js';
import {createSinkModel, pchip, PolarEvaluationError} from '../../assets/js/flight/polar.js';
import {simulateLeg, sinkRate, airGround, GRAVITY} from '../../assets/js/flight/core.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {checkpointRows} from '../../assets/js/ui/checkpoints.js';
import {altitudeTimeSeries, chartSvg} from '../../assets/js/ui/chart.js';
import {loadProject, parseProject, serializeProject, validateProject} from '../../assets/js/io/project.js';
const near=(a,b,tol=1e-8)=>assert.ok(Math.abs(a-b)<=tol,`${a} != ${b}`);
const performance=(E=34,V=25)=>({model:'quadratic-drag',quadratic:{bestLD:E,bestLDAtMps:V},polarTable:null});
const synthetic=(lo=15,hi=45,count=16)=>{
  const b=createSinkModel(performance());
  return Array.from({length:count},(_,i)=>{
    const speedMps=lo+(hi-lo)*i/(count-1);
    return {speedMps,sinkMps:b.rate(speedMps,0)};
  });
};
const polar=(points=synthetic())=>({model:'polar-table',quadratic:{bestLD:34,bestLDAtMps:25},polarTable:{points}});
const state=(speed=25,heading=140)=>{const g=airGround(heading,speed,{fromDeg:270,speedMps:0});return {e:0,n:0,h:500,headingDeg:heading,trackDeg:g.track,bankDeg:0,speedMps:speed,groundMps:g.speed,tS:0,sM:0};};

test('v0.1.13 retains Menuma and 1200ft reference altitude, with model B and schema 1.2',()=>{
  const p=createDemo();
  assert.equal(p.schemaVersion,'1.2');assert.equal(p.performance.model,'quadratic-drag');
  near(p.performance.quadratic.bestLD,34);near(p.performance.quadratic.bestLDAtMps,25);
  near(p.departure.altitudeMSLM-p.landing.elevationMSLM,1200*.3048);
  near(p.landing.point.latDeg,36.211388889,1e-10);
});

test('B is the exact v0.1.11 quadratic formula at several speeds and banks',()=>{
  const p=performance(), m=createSinkModel(p), E=34, vg=25;
  for (const v of [15,20,25,30,50]) for(const bank of [0,15,30,60]) {
    const expected=v**3/(2*E*vg**2)+vg**2/(2*E*v*Math.cos(bank*Math.PI/180)**2);
    near(m.rate(v,bank),expected,1e-12);
    near(sinkRate(v,bank,p),expected,1e-12);
  }
});

test('PCHIP passes given nodes, keeps shape on monotonic intervals and avoids overshoot',()=>{
  const points=[{speedMps:12,sinkMps:1.4},{speedMps:19,sinkMps:1.0},{speedMps:27,sinkMps:.75},{speedMps:41,sinkMps:.9},{speedMps:50,sinkMps:1.7}];
  const f=pchip(points);
  for(const p of points)near(f(p.speedMps),p.sinkMps);
  for(let i=0;i<points.length-1;i++) for(let j=0;j<=100;j++){
    const val=f(points[i].speedMps+(points[i+1].speedMps-points[i].speedMps)*j/100);
    assert.ok(val>=Math.min(points[i].sinkMps,points[i+1].sinkMps)-1e-12);
    assert.ok(val<=Math.max(points[i].sinkMps,points[i+1].sinkMps)+1e-12);
  }
});

test('C accepts exact table edges, rejects extrapolation, without substituting model B',()=>{
  const points=synthetic(15,45,10), m=createSinkModel(polar(points));
  near(m.rate(15,0),points[0].sinkMps);
  near(m.rate(45,0),points.at(-1).sinkMps);
  assert.throws(()=>m.rate(14.9,0),e=>e instanceof PolarEvaluationError&&e.code==='POLAR_SPEED_OUT_OF_RANGE');
  assert.throws(()=>m.rate(45.1,0),e=>e.code==='POLAR_SPEED_OUT_OF_RANGE');
});

test('C uses bank speed V/sqrt(n), and a synthetic B polar approximates the B analytic bank expression',()=>{
  const m=createSinkModel(polar(synthetic(12,60,49)));
  const b=createSinkModel(performance());
  near(m.rate(25,0),b.rate(25,0),1e-12);
  for(const bank of [0,15,30,45,60]){
    const n=1/Math.cos(bank*Math.PI/180);
    const v=28;
    near(m.rate(v,bank),n**1.5*pchip(synthetic(12,60,49))(v/Math.sqrt(n)),1e-10);
    near(m.rate(v,bank),b.rate(v,bank),.002);
  }
});

test('C reports out-of-range at 90km/h / 60 degrees when table starts at 70km/h',()=>{
  const points=synthetic(70/3.6,45,8);
  assert.throws(()=>createSinkModel(polar(points)).rate(25,60),e=>e.code==='POLAR_SPEED_OUT_OF_RANGE');
  const leg={...makeTurn('right',90),bankDeg:60};
  const r=simulateLeg(state(),leg,{speedMps:25,wind:{fromDeg:270,speedMps:0},performance:polar(points)},1,'departure');
  assert.ok(r.notices.some(n=>n.code==='POLAR_SPEED_OUT_OF_RANGE'));
  assert.equal(r.end.altitudeValid,false);
  assert.ok(Number.isFinite(r.end.e) && Number.isFinite(r.end.n));
});

test('C route and Bezier produce error diagnosis without invented checkpoint heights or margin',()=>{
  const p=createDemo();p.performance=polar(synthetic(70/3.6,42,8));
  p.departureLegs=[{...makeTurn('right',90),bankDeg:60}];
  p.arrivalLegs=trafficPreset('left');
  const plan=buildPlan(p);
  assert.equal(plan.valid,false);
  assert.ok(plan.notices.some(n=>n.code==='POLAR_SPEED_OUT_OF_RANGE'));
  assert.equal(plan.departure.end.altitudeValid,false);
  assert.ok(plan.departure.samples.every(s=>Number.isFinite(s.h) && Number.isFinite(s.e)));
  assert.ok(checkpointRows(plan,p).some(r=>r.kind==='height' && r.relativeM===undefined));
  assert.ok(!plan.checkpoints.some(cp=>cp.marginM!==undefined));
  assert.ok(!chartSvg(plan,p).includes('NaN'));
  assert.ok(altitudeTimeSeries(plan).durationS>0);
});

test('C and B settings survive switching and JSON roundtrip',()=>{
  const p=createDemo();p.performance=polar(synthetic());
  p.performance.quadratic.bestLD=37;
  const saved=parseProject(serializeProject(p));
  assert.deepEqual(saved,p);
  saved.performance.model='quadratic-drag';
  assert.equal(saved.performance.quadratic.bestLD,37);
  assert.deepEqual(saved.performance.polarTable,p.performance.polarTable);
  assert.doesNotThrow(()=>validateProject(saved));
});

test('C input rejects missing, repeated, unordered and non-positive points',()=>{
  const p=createDemo();p.performance=polar(synthetic(15,45,4));
  const broken=structuredClone(p);broken.performance.polarTable=null;
  assert.throws(()=>validateProject(broken),/POLAR_DATA_MISSING/);
  broken.performance.polarTable={points:[{speedMps:10,sinkMps:1}]};
  assert.throws(()=>validateProject(broken),/4〜64/);
  broken.performance.polarTable={points:structuredClone(p.performance.polarTable.points)};
  broken.performance.polarTable.points[1].speedMps=broken.performance.polarTable.points[0].speedMps;
  assert.throws(()=>validateProject(broken),/重複のない昇順/);
  broken.performance.polarTable.points=structuredClone(p.performance.polarTable.points);
  broken.performance.polarTable.points[2].sinkMps=0;
  assert.throws(()=>validateProject(broken),/正の有限値/);
});

function oldModelB(){
  const p=createDemo();p.schemaVersion='1.0';p.engineVersion='0.1.11-test';
  p.performance={model:'quadratic-drag',fixedLD:28,bestLD:31,bestLDAtMps:22};
  p.defaultSpeedMps=24.5;
  return p;
}

test('unsupported old model is rejected without automatic conversion or modifying input',()=>{
  const old=oldModelB();old.performance.model='fixed-ld';
  const before=structuredClone(old);
  assert.throws(()=>loadProject(old),/簡易抗力モデルBのみ/);
  assert.throws(()=>parseProject(JSON.stringify(old)),/簡易抗力モデルBのみ/);
  assert.deepEqual(old,before);
});

test('legacy B is migrated automatically with identical analytic sink and unchanged positions/legs',()=>{
  const old=oldModelB();
  old.departureLegs[0].speedOverrideMps=30;
  const modern=loadProject(old);
  near(createSinkModel(modern.performance).rate(27,30),27**3/(2*31*22**2)+22**2/(2*31*27*Math.cos(Math.PI/6)**2));
  assert.deepEqual(modern.departure,old.departure);
  assert.deepEqual(modern.arrivalLegs,old.arrivalLegs);
  assert.equal(modern.departureLegs[0].speedOverrideMps,30);
  assert.equal(parseProject(serializeProject(modern)).schemaVersion,'1.2');
});

test('old fixed-L/D model cannot be saved as the current 1.2 format',()=>{
  const p=createDemo();p.performance.model='fixed-ld';
  assert.throws(()=>serializeProject(p),/計算モデル/);
});

test('flight UI and PDF only display model B/C, not a default airfield/aircraft name',()=>{
  const main=readFileSync(new URL('../src/main.ts',import.meta.url),'utf8');
  const print=readFileSync(new URL('../src/ui/print.ts',import.meta.url),'utf8');
  assert.match(main,/C：速度別ポーラーモデル/);
  assert.match(print,/C：速度別ポーラー/);
  assert.doesNotMatch(main,/A：固定L\/D|妻沼滑空場|ASK 21/);
  assert.doesNotMatch(print,/A：固定L\/D|妻沼滑空場|ASK 21/);
});
