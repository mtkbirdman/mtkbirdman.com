import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo,makeStraight,makeTurn,trafficPreset,ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {effectiveLegSpeedMps,effectiveStraightDirectionDeg} from '../../assets/js/ui/leg-display.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {serializeProject,parseProject} from '../../assets/js/io/project.js';
const read=file=>readFileSync(new URL(file,import.meta.url),'utf8');
const close=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<tol,`${a} differs from ${b}`);

test('v0.1.10 versions agree across source, WordPress, print, and documentation',()=>{
 assert.equal(ENGINE_VERSION,'0.1.20-test');
 assert.match(read('../../glider-landing-planner.php'),/Version:\s*0\.1\.20/);
 assert.match(read('../src/main.ts'),/テスト版 v0\.1\.20/);
 assert.match(read('../src/ui/print.ts'),/ver=0\.1\.20/);
 assert.match(read('../../README.md'),/v0\.1\.20/);
});

test('v0.1.9 displays inherited baseline speed without converting it to an override',()=>{
 const p=createDemo();const leg=makeStraight(300);p.departureLegs=[leg];
 close(effectiveLegSpeedMps(p,leg),90/3.6);
 assert.equal(leg.speedOverrideMps,null);
 p.defaultSpeedMps=105/3.6;
 close(effectiveLegSpeedMps(p,leg),105/3.6);
 assert.equal(parseProject(serializeProject(p)).departureLegs[0].speedOverrideMps,null);
 leg.speedOverrideMps=90/3.6;
 p.defaultSpeedMps=120/3.6;
 close(effectiveLegSpeedMps(p,leg),90/3.6);
});

test('v0.1.9 forward inherited target shows actual heading after preceding turn',()=>{
 const p=createDemo();p.departure.headingDeg=272;p.wind={fromDeg:0,speedMps:5};
 p.departureLegs=[makeStraight(150),makeTurn('right',90),makeStraight(150)].map(l=>({...l,guidance:'heading'}));
 const plan=buildPlan(p);
 const first=effectiveStraightDirectionDeg(p,plan,'departure',p.departureLegs[0]);
 const later=effectiveStraightDirectionDeg(p,plan,'departure',p.departureLegs[2]);
 close(first,272,0.01);close(later,2,0.01);
 assert.equal(p.departureLegs[2].directionDeg,null);
});

test('v0.1.9 ground-track target shows track, not heading, with crosswind',()=>{
 const p=createDemo();p.wind={fromDeg:270,speedMps:4};p.departure.headingDeg=0;
 const leg={...makeStraight(100),guidance:'ground-track'};p.departureLegs=[leg];
 const plan=buildPlan(p);
 const seg=plan.departure.segments[0];
 close(effectiveStraightDirectionDeg(p,plan,'departure',leg),seg.samples[0].trackDeg);
 assert.notEqual(seg.samples[0].headingDeg,seg.samples[0].trackDeg);
});

test('v0.1.9 arrival leg inherits from actual reverse-time computed state',()=>{
 const p=createDemo();p.wind={fromDeg:0,speedMps:4};
 p.arrivalLegs=[makeStraight(250),makeTurn('right',90),makeStraight(200)].map(l=>({...l,guidance:'heading'}));
 const plan=buildPlan(p);
 for(const leg of p.arrivalLegs.filter(l=>l.type==='straight')){
   const seg=plan.arrival.segments.find(s=>s.leg.id===leg.id);
   close(effectiveStraightDirectionDeg(p,plan,'arrival',leg),seg.samples[0].headingDeg);
 }
 const last=p.arrivalLegs.at(-1);
 last.directionDeg=23;
 close(effectiveStraightDirectionDeg(p,buildPlan(p),'arrival',last),23);
});

test('presets and manual additions default to ground-track',()=>{
 for(const side of ['left','right']){
   const legs=trafficPreset(side);
   assert.equal(legs.length,4);
   assert.ok(legs.every(l=>l.guidance==='ground-track'));
   assert.deepEqual(legs.map(l=>l.checkpointLabel),['ダウンウィンド開始','ベースターン開始','ベースレグ開始','ファイナルターン開始']);
 }
 assert.equal(makeStraight().guidance,'ground-track');
 assert.equal(makeTurn().guidance,'ground-track');
});

test('v0.1.9 import of legacy heading-mode preset does not silently rewrite it',()=>{
 const old=createDemo();old.engineVersion='0.1.8-test';
 old.arrivalLegs=trafficPreset('right').map(l=>({...l,guidance:'heading'}));
 const copy=parseProject(serializeProject(old));
 assert.ok(copy.arrivalLegs.every(l=>l.guidance==='heading'));
 assert.equal(copy.engineVersion,'0.1.8-test');
});

test('v0.1.9 form and report labels use updated wording without changing data keys',()=>{
 const main=read('../src/main.ts');const print=read('../src/ui/print.ts');
 for(const word of ['基準対気速度','飛行方式','進行方向','滑空速度','旋回速度'])assert.ok(main.includes(word),word);
 for(const word of ['共通対気速度','風補正方式','指定方位（空欄＝継承）','多周旋回可','速度を空欄にした場合'])assert.ok(!main.includes(word),word);
 assert.match(main,/effectiveLegSpeedMps\(this\.project,\s*leg\)/);
 assert.match(main,/effectiveStraightDirectionDeg\(this\.project,\s*this\.plan,\s*side,\s*leg\)/);
 assert.match(main,/this\.control\(base\s*\+\s*'\.speedOverrideMps'[^\n]*speedInForce\)/);
 assert.match(main,/this\.control\(base\s*\+\s*'\.directionDeg'[^\n]*directionInForce\)/);
 assert.match(print,/<th>飛行方式<\/th>/);
 assert.match(print,/基準対気速度/);
});
