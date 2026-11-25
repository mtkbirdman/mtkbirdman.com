import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo, makeTurn} from '../../assets/js/model/defaults.js';
import {airGround, closestTurnAngle, simulateLeg} from '../../assets/js/flight/core.js';
import {buildPlan, rememberResolvedTurnAngles} from '../../assets/js/planning/plan.js';
import {effectiveTurnValues} from '../../assets/js/ui/leg-display.js';
import {loadProject, validateProject, parseProject, serializeProject} from '../../assets/js/io/project.js';
import {reportHtml} from '../../assets/js/ui/print.js';
import {readFileSync} from 'node:fs';
const close=(a,b,tol=1e-5)=>assert.ok(Math.abs(a-b)<tol,`expected ${a} ≈ ${b}`);
const makeTarget=(side='right',reference=90,target=180,guidance='heading')=>({...makeTurn(side,reference),angleInputMode:'target',targetDirectionDeg:target,guidance});
const wind={fromDeg:270,speedMps:5};
const state=(heading=90,speedMps=25,weather={fromDeg:270,speedMps:0})=>{const g=airGround(heading,speedMps,weather);return {e:0,n:0,h:500,headingDeg:heading,trackDeg:g.track,bankDeg:0,speedMps,groundMps:g.speed,tS:0,sM:0};};
const context=(weather={fromDeg:270,speedMps:0})=>({speedMps:25,wind:weather,performance:createDemo().performance});

test('v0.1.13 starts with angle input and schema 1.3 with the new reference altitude',()=>{
  const p=createDemo(),turn=makeTurn();
  assert.equal(p.schemaVersion,'1.3');assert.equal(turn.angleInputMode,'angle');assert.equal(turn.targetDirectionDeg,null);
  close(p.departure.altitudeMSLM-p.landing.elevationMSLM,1200*.3048,1e-10);
  assert.deepEqual(p.performance,{bestLD:34,bestLDAtMps:25});
});

test('right and left target conversion is correct for both planning directions',()=>{
  close(closestTurnAngle(180,90,'right',false,90),90);
  close(closestTurnAngle(0,90,'left',false,90),90);
  close(closestTurnAngle(0,90,'right',true,90),90);
  close(closestTurnAngle(180,90,'left',true,90),90);
});

test('multi-revolution target solution stays nearest positive and ties choose the smaller angle',()=>{
  close(closestTurnAngle(180,90,'right',false,450),450);
  close(closestTurnAngle(180,90,'right',false,810),810);
  close(closestTurnAngle(180,90,'right',false,270),90); // tie: 90 and 450
  close(closestTurnAngle(270,90,'right',false,90),180);
  close(closestTurnAngle(90,90,'right',false,360),360);
  close(closestTurnAngle(90,90,'right',false,720),720);
  close(closestTurnAngle(90,90,'right',true,360),360);
});

test('heading-driven departure target resolves 450 degrees without collapsing to 90',()=>{
  const leg=makeTarget('right',450,180);
  const r=simulateLeg(state(),leg,context(),1,'departure');
  close(r.effectiveLeg.angleDeg,450);
  close(r.end.headingDeg,180);
  assert.equal(leg.angleDeg,450); // caller input stays untouched
  close(r.samples[0].headingDeg,90);
});

test('arrival target constrains the pre-turn heading while keeping real-flight turn direction',()=>{
  const leg=makeTarget('right',450,0);
  const r=simulateLeg(state(),leg,context(),-1,'arrival');
  close(r.effectiveLeg.angleDeg,450);
  close(r.end.headingDeg,0);
  close(r.samples[0].headingDeg,90);
});

test('departure target follows upstream heading changes; angle-driven turn keeps its angle',()=>{
  const p=createDemo();p.departure.headingDeg=90;
  const leg=makeTarget('right',90,180);
  p.departureLegs=[leg];
  close(buildPlan(p).departure.segments[0].leg.angleDeg,90);
  p.departure.headingDeg=100;
  close(buildPlan(p).departure.segments[0].leg.angleDeg,80);
  leg.angleInputMode='angle';leg.targetDirectionDeg=null;leg.angleDeg=90;
  const plan=buildPlan(p);
  close(plan.departure.segments[0].leg.angleDeg,90);
  close(effectiveTurnValues(plan,'departure',leg).boundaryDeg,190);
});

test('ground-track target uses track at leg override speed after wind correction',()=>{
  const initial=state(0,20,wind);const ctx=context(wind);
  const oldTrack=initial.trackDeg;
  const turn={...makeTarget('right',90,135,'ground-track'),speedOverrideMps:30};
  const newTrack=airGround(0,30,wind).track;
  const r=simulateLeg(initial,turn,ctx,1,'departure');
  close(r.effectiveLeg.angleDeg,closestTurnAngle(135,newTrack,'right',false,90));
  assert.ok(Math.abs(closestTurnAngle(135,oldTrack,'right',false,90)-r.effectiveLeg.angleDeg)>1);
  close(r.end.trackDeg,135,0.025);
  assert.ok(Math.abs(r.end.headingDeg-90)>1);
});

test('arrival ground-track target is resolved from final-side boundary after speed correction',()=>{
  const initial=state(140,20,wind);const ctx=context(wind);
  const leg={...makeTarget('right',90,45,'ground-track'),speedOverrideMps:30};
  const actualBoundary=airGround(140,30,wind).track;
  const r=simulateLeg(initial,leg,ctx,-1,'arrival');
  close(r.effectiveLeg.angleDeg,closestTurnAngle(45,actualBoundary,'right',true,90));
  close(r.end.trackDeg,45,.025);
});

test('ground-track target does not bypass unreachable-turn solver diagnostics',()=>{
  const strongWind={fromDeg:270,speedMps:30};
  const r=simulateLeg(state(0,25,strongWind),makeTarget('right',90,0,'ground-track'),context(strongWind),1,'departure');
  assert.ok(r.notices.some(n=>n.code==='TURN_TARGET_NOT_REACHED'));
});

test('reverse-resolved arrival angle is identical to forward predicted arrival angle',()=>{
  const p=previouslyPlannedProject();p.arrivalLegs=[makeTarget('right',450,50)];
  const plan=buildPlan(p);
  close(plan.arrival.segments[0].leg.angleDeg,450);
  close(plan.arrival.start.headingDeg,50);
  const forwardSamples=plan.forwardAll.filter(s=>s.segmentId===p.arrivalLegs[0].id);
  assert.ok(forwardSamples.length>100);
  const last=forwardSamples.at(-1);
  close(last.e,plan.arrival.end.e,1e-5);
  close(last.n,plan.arrival.end.n,1e-5);
  close(last.headingDeg,plan.arrival.end.headingDeg,1e-5);
});

test('target-mode stored anchor is updated from solved segments without touching target bearing',()=>{
  const p=createDemo();p.departure.headingDeg=100;p.departureLegs=[makeTarget('right',90,180)];
  const plan=buildPlan(p);
  rememberResolvedTurnAngles(p,plan);
  assert.equal(p.departureLegs[0].angleDeg,80);
  assert.equal(p.departureLegs[0].targetDirectionDeg,180);
  p.departureLegs[0].targetDirectionDeg=190;
  const next=buildPlan(p);
  close(next.departure.segments[0].leg.angleDeg,90);
});


test('schema 1.3 target mode roundtrips; malformed target is rejected, not normalized',()=>{
  const p=createDemo();p.departureLegs=[makeTarget('right',810,180)];
  const parsed=parseProject(serializeProject(p));
  assert.deepEqual(parsed,p);
  for(const bad of [null,-1,361,'90',NaN]){
    const broken=structuredClone(p);broken.departureLegs[0].targetDirectionDeg=bad;
    assert.throws(()=>validateProject(broken),/目標方位/);
  }
});


test('PDF prints the resolved turn angle, without altering the report layout',()=>{
  globalThis.document={querySelector:()=>null,querySelectorAll:()=>[]};
  const p=createDemo();p.departure.headingDeg=100;p.departureLegs=[makeTarget('right',90,180)];
  const plan=buildPlan(p);
  const mapCopy={querySelectorAll:()=>[],removeAttribute(){},classList:{add(){}},style:{},outerHTML:'<div class="leaflet"></div>'};
  const mapEl={clientWidth:900,clientHeight:420,cloneNode(){return mapCopy}};
  const html=reportHtml({project:p,plan,mapEl,metricsHtml:'',diagnosticsHtml:'',title:'GLP',engineVersion:'0.1.17-test'});
  assert.match(html,/右旋回 80°、バンク 30°/);
  assert.doesNotMatch(html,/右旋回 90°、バンク 30°/);
  assert.match(html,/飛行経路図/);assert.match(html,/高度プロファイル/);
});

globalThis.document={querySelector:()=>null,querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');
globalThis.HTMLInputElement=class {};
globalThis.HTMLSelectElement=class {constructor(path,value){this.dataset={path};this.value=value}};
function mockedApp(){
  const app=Object.create(App.prototype);app.project=createDemo();app.plan=buildPlan(app.project);
  app.history=[];app.future=[];app.editRevision=0;
  app.selected=null;app.notice='';app.status='';app.render=()=>{};app.renderResults=()=>{};
  app.renderSidebar=()=>{};app.scheduleSave=()=>{};
  return app;
}

test('PC and mobile shared editor displays both fields with only selected one writable',()=>{
  const app=mockedApp();app.project.departure.headingDeg=90;
  const leg=makeTarget('right',450,180);app.project.departureLegs=[leg];app.plan=buildPlan(app.project);
  const html=app.legEditor('departure',leg,0);
  assert.match(html,/目標方位を指定/);assert.match(html,/旋回後の目標方位/);
  assert.doesNotMatch(html,/（自動計算）/);
  assert.match(html,/readonly value="450"/);
  assert.match(html,/data-path="departureLegs\.0\.targetDirectionDeg"/);
  const angleLeg=makeTurn('left',270);app.project.arrivalLegs=[angleLeg];app.plan=buildPlan(app.project);
  const arrivalHtml=app.legEditor('arrival',angleLeg,0);
  assert.match(arrivalHtml,/旋回前の目標方位/);
  assert.match(arrivalHtml,/data-path="arrivalLegs\.0\.angleDeg"/);
  assert.doesNotMatch(arrivalHtml,/（自動計算）/);
});

test('changing modes initializes the target from actual track and switching back freezes the effective angle',()=>{
  const app=mockedApp();app.project.departure.headingDeg=90;
  app.project.departureLegs=[makeTurn('right',450)];app.plan=buildPlan(app.project);
  app.handleChange(new HTMLSelectElement('departureLegs.0.angleInputMode','target'));
  assert.equal(app.project.departureLegs[0].angleInputMode,'target');
  close(app.project.departureLegs[0].targetDirectionDeg,180);
  close(app.project.departureLegs[0].angleDeg,450);
  app.handleChange(new HTMLSelectElement('departureLegs.0.angleInputMode','angle'));
  assert.equal(app.project.departureLegs[0].angleInputMode,'angle');
  assert.equal(app.project.departureLegs[0].targetDirectionDeg,null);
  close(app.project.departureLegs[0].angleDeg,450);
  assert.equal(app.history.length,2);
});

test('editing target bearing updates effective angle and undo/redo restores mode and anchor',()=>{
  const app=mockedApp();app.project.departure.headingDeg=90;
  app.project.departureLegs=[makeTarget('right',450,180)];app.plan=buildPlan(app.project);
  // mimic numeric field input change without DOM
  assert.equal(app.edit(p=>{p.departureLegs[0].targetDirectionDeg=200}),true);
  close(app.project.departureLegs[0].angleDeg,470);
  assert.equal(app.project.departureLegs[0].targetDirectionDeg,200);
  app.handleAction('undo',{});
  close(app.project.departureLegs[0].angleDeg,450);
  assert.equal(app.project.departureLegs[0].targetDirectionDeg,180);
  app.handleAction('redo',{});
  close(app.project.departureLegs[0].angleDeg,470);
});

test('optimization signature includes mode and target, and editing recomputes the connection',()=>{
  const app=mockedApp();app.project.departureLegs=[makeTurn('right',90)];app.plan=buildPlan(app.project);
  const before=app.optimizationSignature(app.project);
  app.edit(p=>{p.departureLegs[0].angleInputMode='target';p.departureLegs[0].targetDirectionDeg=180;});
  assert.notEqual(app.optimizationSignature(app.project),before);
  const old=app.plan.departure.segments[0].leg.angleDeg;
  app.edit(p=>{p.departureLegs[0].targetDirectionDeg=200;});
  assert.notEqual(app.plan.departure.segments[0].leg.angleDeg,old);
});

test('single sink model retains the WordPress shortcode and magnetic convention',()=>{
  const root=new URL('../',import.meta.url);
  const php=readFileSync(new URL('../glider-landing-planner.php',root),'utf8');
  assert.match(php,/Version:\s*0\.1\.17/);assert.match(php,/glider_landing_planner/);
  const p=createDemo();assert.equal(p.display.magneticDeclinationDeg,0);
  assert.deepEqual(p.performance,{bestLD:34,bestLDAtMps:25});
});
