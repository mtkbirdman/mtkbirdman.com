import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo,makeTurn,ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {parseProject,serializeProject} from '../../assets/js/io/project.js';
import {simulateLeg,simulatePath,airGround,GRAVITY,TURN_STEP_BUDGET} from '../../assets/js/flight/core.js';
import {magneticToTrue,rad} from '../../assets/js/geometry/geo.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {chooseConnector} from '../../assets/js/planning/connector.js';
import {altitudeTimeSeries} from '../../assets/js/ui/chart.js';
const read=(p)=>readFileSync(new URL(p,import.meta.url),'utf8');
const wind={fromDeg:270,speedMps:3};
const start=(w={fromDeg:0,speedMps:0})=>{const a=airGround(0,20,w);return {e:0,n:0,h:1000,headingDeg:0,trackDeg:a.track,bankDeg:0,speedMps:20,groundMps:a.speed,tS:0,sM:0}};
const ctx=(w={fromDeg:0,speedMps:0})=>({speedMps:20,wind:w,performance:{model:'fixed-ld',fixedLD:30,bestLD:30,bestLDAtMps:20}});
const near=(v,x,tol)=>assert.ok(Math.abs(v-x)<tol,`${v} differs from ${x}`);
test('v0.1.6 positive unbounded turn angles survive JSON with multi-turn records',()=>{
 const p=createDemo();assert.equal(ENGINE_VERSION,'0.1.8-test');
 p.departureLegs=[makeTurn('right',1080),makeTurn('left',1800),makeTurn('right',1000000)];
 const q=parseProject(serializeProject(p));assert.deepEqual(q.departureLegs.map(x=>x.angleDeg),[1080,1800,1000000]);
 p.departureLegs[0].angleDeg=0;assert.throws(()=>parseProject(JSON.stringify(p)),/旋回角/);
});
test('v0.1.6 five complete still-air heading circles close and integrate time/loss',()=>{
 const l=makeTurn('right',1800),r=simulateLeg(start(),l,ctx(),1,'departure');
 assert.equal(r.notices.some(x=>x.level==='error'),false);
 near(r.end.e,0,.03);near(r.end.n,0,.03);near(r.end.headingDeg,0,.0001);
 const t=rad(1800)*20/(GRAVITY*Math.tan(rad(30)));
 near(r.end.tS,t,.01);near(1000-r.end.h,(20/30)*t,.01);
 assert.ok(r.samples.length>6000);
});
test('v0.1.6 three complete ground-track turns reach cumulative 1080 with wind',()=>{
 const l={...makeTurn('right',1080),guidance:'ground-track'},r=simulateLeg(start(wind),l,ctx(wind),1,'departure');
 assert.equal(r.notices.some(x=>x.level==='error'),false,JSON.stringify(r.notices));
 near(r.end.trackDeg,r.samples[0].trackDeg,0.03);
 assert.ok(r.end.sM>0&&r.end.tS>0);
});
test('v0.1.6 extreme angles terminate safely without looping indefinitely',()=>{
 const angle=1000000,l=makeTurn('right',angle),r=simulateLeg(start(),l,ctx(),1,'departure');
 assert.ok(TURN_STEP_BUDGET>1000);
 assert.ok(r.notices.some(x=>x.code==='TURN_COMPLEXITY_LIMIT'&&x.level==='error'));
 assert.ok(r.samples.length<=2);
});
test('v0.1.6 magnetic declination changes only overlay true orientation',()=>{
 near(magneticToTrue(0,7),7,1e-10);near(magneticToTrue(0,-7),353,1e-10);
 near(magneticToTrue(350,20),10,1e-10);
 const p=createDemo(),before=buildPlan(p);p.display.magneticDeclinationDeg=-7;
 const after=buildPlan(p);near(before.departure.end.e,after.departure.end.e,1e-9);near(before.departure.end.n,after.departure.end.n,1e-9);
});
test('v0.1.6 old JSON migrates shared opacity + intervals + connector mode',()=>{
 const p=createDemo();delete p.display.referenceLineOpacity;
 delete p.display.distanceRingIntervalM;delete p.display.bearingLineIntervalDeg;delete p.display.magneticDeclinationDeg;
 p.display.distanceRingOpacity=0.3;p.display.bearingLineOpacity=0.8;
 delete p.connector.adjustmentMode;
 const q=parseProject(JSON.stringify(p));near(q.display.referenceLineOpacity,0.3,1e-9);
 assert.equal(q.display.distanceRingIntervalM,250);assert.equal(q.display.bearingLineIntervalDeg,45);assert.equal(q.display.magneticDeclinationDeg,0);
 assert.equal(q.connector.adjustmentMode,'auto');
 assert.deepEqual(parseProject(serializeProject(q)).display.referenceLineOpacity,q.display.referenceLineOpacity);
});
test('v0.1.6 range checks reject invalid ring spacing/compass interval/declination/opacity',()=>{
 const p=createDemo();p.display.distanceRingIntervalM=2;assert.throws(()=>parseProject(JSON.stringify(p)),/距離リング間隔/);
 p.display.distanceRingIntervalM=250;p.display.bearingLineIntervalDeg=25;assert.throws(()=>parseProject(JSON.stringify(p)),/方位線間隔/);
 p.display.bearingLineIntervalDeg=45;p.display.magneticDeclinationDeg=90;assert.throws(()=>parseProject(JSON.stringify(p)),/磁気偏角/);
 p.display.magneticDeclinationDeg=0;p.display.referenceLineOpacity=-.1;assert.throws(()=>parseProject(JSON.stringify(p)),/不透明度/);
});
test('v0.1.6 editable auto/manual Bezier mode is saved and valid',()=>{
 const p=createDemo();assert.equal(p.connector.adjustmentMode,'auto');
 p.connector.adjustmentMode='manual';assert.equal(parseProject(serializeProject(p)).connector.adjustmentMode,'manual');
 p.connector.adjustmentMode='invalid';assert.throws(()=>parseProject(JSON.stringify(p)),/Bezier調整モード/);
});
test('v0.1.6 auto-adjust is triggered only on flight geometry/conditions, manual handle changes disable auto',()=>{
 const code=read('../src/main.ts');
 assert.match(code,/private optimizationSignature\(/);
 assert.match(code,/modeChanged=previous\.connector\?\.adjustmentMode!=='auto'/);
 assert.match(code,/shapeChanged=this\.optimizationSignature\(previous\)!==this\.optimizationSignature\(p\)/);
 assert.match(code,/chooseConnector\(preliminary\.departure\.end,preliminary\.arrivalStart,p\)/);
 assert.match(code,/p\.connector\.adjustmentMode='manual'/);
 assert.doesNotMatch(code,/this\.control\('connector\.handleStartM'/);
 assert.match(code,/接続曲線の調整/);
});
test('v0.1.6 Bezier search reoptimizes after a route geometry edit',()=>{
 const p=createDemo(),r=buildPlan(p),a=chooseConnector(r.departure.end,r.arrivalStart,p);
 p.departureLegs[0].lengthM+=300;
 const r2=buildPlan(p),b=chooseConnector(r2.departure.end,r2.arrivalStart,p);
 assert.ok(Number.isFinite(a.handleStartM)&&Number.isFinite(b.handleEndM));
 assert.ok(a.handleStartM!==b.handleStartM||a.handleEndM!==b.handleEndM,'geometry edit should modify candidate handles');
});
test('v0.1.6 results include flight time and margin; diagnostic contains bank, PDF heading not duplicated',()=>{
 const main=read('../src/main.ts'),print=read('../src/ui/print.ts'),css=read('../../assets/css/print.css');
 assert.match(main,/推定飛行時間/);assert.match(main,/高度余裕<\/small>/);
 assert.match(main,/Bezier接続最大必要バンク/);
 assert.match(main,/diag\.querySelector\(sel\)\?\.outerHTML/);
 assert.match(print,/<h2>計算上の診断<\/h2>/);
 assert.doesNotMatch(print,/diagnosticsHtml\.innerHTML/);
 assert.match(css,/glp-report-map-frame\{[^\n]*background:#fff/);
 assert.match(print,/width:\$\{Math\.ceil\(width\*scale\)\}px/);
});
test('v0.1.6 ring count and magnetic bearing options are in UI and map source',()=>{
 const main=read('../src/main.ts'),map=read('../src/ui/map.ts');
 assert.match(main,/display\.referenceLineOpacity/);assert.match(main,/display\.distanceRingIntervalM/);
 assert.match(main,/display\.bearingLineIntervalDeg/);assert.match(main,/display\.magneticDeclinationDeg/);
 assert.match(map,/magneticToTrue\(magnetic,project\.display\.magneticDeclinationDeg\)/);
 assert.match(map,/radius\+=step/);
 assert.match(map,/magnetic\+=project\.display\.bearingLineIntervalDeg/);
});
