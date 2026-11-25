import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo,trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {loadProject,validateProject,parseProject,serializeProject} from '../../assets/js/io/project.js';
import {escapeHtml} from '../../assets/js/ui/escape-html.js';
import {fmtAlt} from '../../assets/js/ui/chart.js';
import {heightForDisplay} from '../../assets/js/ui/height-units.js';

// Import the browser entry point with no app roots; no DOM library or production
// test hooks are needed. Its TypeScript-private methods remain normal JS methods.
globalThis.document={querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');


function bareEditor(){
 const app=Object.create(App.prototype);
 app.project=previouslyPlannedProject();
 app.plan=buildPlan(app.project);
 app.history=[];
 app.future=[];
 app.editRevision=0;
 app.notice='';
 app.status='';
 app.renders=0;
 app.saves=0;
 app.render=()=>app.renders++;
 app.scheduleSave=()=>app.saves++;
 return app;
}

test('v0.1.10 display-only and name changes retain the exact calculated plan',()=>{
 const app=bareEditor(),plan=app.plan;
 assert.equal(app.edit(p=>p.display.altitudeUnit='m'),true);
 assert.equal(app.plan,plan,'display changes must not re-run the flight calculation');
 assert.equal(app.edit(p=>p.name='テスト計画'),true);
 assert.equal(app.plan,plan,'renaming a project must not re-run the flight calculation');
 assert.equal(app.history.length,2);
 assert.equal(app.saves,2);
});

test('v0.1.10 checkpoint rename updates labels without rerunning Bezier optimization',()=>{
 const app=bareEditor();
 const before=structuredClone(app.project.connector);
 assert.equal(app.edit(p=>p.departureLegs[0].checkpointLabel='独立したラベル'),true);
 assert.equal(app.plan.checkpoints.find(cp=>cp.group==='departure'&&cp.legId)?.label,'独立したラベル');
 assert.deepEqual(app.project.connector,before);
});

test('v0.1.10 geometry edits still rebuild and reoptimize an auto connector',()=>{
 const app=bareEditor(),oldPlan=app.plan,oldLength=app.project.departureLegs[0].lengthM;
 assert.equal(app.edit(p=>p.departureLegs[0].lengthM+=300),true);
 assert.notEqual(app.plan,oldPlan);
 assert.equal(app.project.departureLegs[0].lengthM,oldLength+300);
 assert.notDeepEqual(app.plan.departure.end,oldPlan.departure.end);
 assert.equal(app.plan.connector.handleStartM,app.project.connector.handleStartM);
});

test('v0.1.10 no-op edit does not populate history or schedule autosave',()=>{
 const app=bareEditor();
 assert.equal(app.edit(p=>{p.display.showHandles=p.display.showHandles;}),false);
 assert.equal(app.history.length,0);
 assert.equal(app.saves,0);
});

test('v0.1.10 undo/redo share the same commit path and keep the exact geometry',()=>{
 const app=bareEditor(),original=app.project.departureLegs[0].lengthM;
 app.edit(p=>p.departureLegs[0].lengthM=original+100);
 const modified=app.plan.departure.end.e;
 app.handleAction('undo',{});
 assert.equal(app.project.departureLegs[0].lengthM,original);
 assert.equal(app.future.length,1);
 app.handleAction('redo',{});
 assert.equal(app.project.departureLegs[0].lengthM,original+100);
 assert.equal(app.plan.departure.end.e,modified);
 assert.equal(app.saves,3);
});

test('v0.1.10 validation has no mutation side effects, including on rejected records',()=>{
 const current=createDemo(),before=structuredClone(current);
 assert.equal(validateProject(current),current);
 assert.deepEqual(current,before);
 const broken=structuredClone(current);
 broken.wind.speedMps=-1;
 assert.throws(()=>validateProject(broken),/風向・風速/);
 assert.deepEqual(broken.wind,{fromDeg:270,speedMps:-1});
 assert.deepEqual(JSON.parse(serializeProject(current)),current);
});

test('v0.1.10 legacy loading creates a migrated copy, preserving the source object',()=>{
 const legacy=createDemo();
 legacy.engineVersion='0.1.3-test';
 legacy.schemaVersion='1.0';
 legacy.performance={model:'quadratic-drag',fixedLD:30,bestLD:30,bestLDAtMps:20};
 legacy.arrivalLegs=trafficPreset('right');
 legacy.arrivalLegs[1].bankDeg=12;
 delete legacy.arrivalLegs[0].checkpointKind;
 delete legacy.display.bearingLineIntervalDeg;
 delete legacy.display.referenceLineOpacity;
 legacy.display.distanceRingOpacity=0.4;
 const original=structuredClone(legacy);
 const result=loadProject(legacy);
 assert.deepEqual(legacy,original,'migration must not mutate the input');
 assert.equal(result.engineVersion,'0.1.16-test');
 assert.equal(result.arrivalLegs[1].bankDeg,30);
 assert.equal(result.arrivalLegs[0].checkpointKind,'downwind-start');
 assert.equal(result.display.bearingLineIntervalDeg,45);
 assert.equal(result.display.referenceLineOpacity,0.4);
 assert.equal(parseProject(JSON.stringify(legacy)).display.referenceLineOpacity,0.4);
});

test('v0.1.10 legacy migration never disguises an invalid checkpoint kind',()=>{
 const bad=createDemo();bad.arrivalLegs=trafficPreset('left');
 bad.arrivalLegs[0].checkpointKind='invalid';
 assert.throws(()=>loadProject(bad),/チェックポイント種別/);
});

test('v0.1.10 one HTML escaping implementation protects project and checkpoint text',()=>{
 assert.equal(escapeHtml(`<script src='x'>&"</script>`),'&lt;script src=&#39;x&#39;&gt;&amp;&quot;&lt;/script&gt;');
 assert.equal(escapeHtml(null),'');
});

test('v0.1.10 height formatting uses the canonical unit conversion',()=>{
 assert.equal(fmtAlt(304.8,'ft'),'1,000 ft');
 assert.ok(Math.abs(heightForDisplay(304.8,'ft')-1000)<1e-9);
 assert.equal(fmtAlt(304.8,'m'),'305 m');
});

test('v0.1.10 final path-angle diagnostic reflects a modified project',()=>{
 const project=createDemo();project.landing.final.pathAngleDeg=4.5;
 const before=buildPlan(project);
 assert.match(before.notices.find(n=>n.code==='FINAL_BRAKE_UNKNOWN').message,/4\.5°/);
 assert.doesNotMatch(before.notices.find(n=>n.code==='FINAL_BRAKE_UNKNOWN').message,/6°/);
});

test('v0.1.10 adding a leg selects it using existing selection handling',()=>{
 const app=bareEditor();
 app.selected=null;
 app.map={redraw:()=>{},fit:()=>{}};
 app.root={querySelectorAll:()=>[],querySelector:()=>null};
 app.renderSidebar=()=>{};
 app.updateCheckpointSelection=()=>{};
 app.openMobilePanel=()=>{};
 const initialCount=app.project.departureLegs.length;
 app.handleAction('add-straight',{side:'departure'});
 assert.equal(app.project.departureLegs.length,initialCount+1);
 assert.equal(app.selected,app.project.departureLegs.at(-1).id);
 assert.equal(app.history.length,1);
});

test('v0.1.10 deleting a selected leg clears the selection and preserves undo',()=>{
 const app=bareEditor();
 const selected=app.project.departureLegs[0].id;
 app.selected=selected;
 app.mobilePanel=null;
 app.map={redraw:()=>{},fit:()=>{}};
 app.root={querySelectorAll:()=>[],querySelector:()=>null};
 app.renderSidebar=()=>{};
 app.updateCheckpointSelection=()=>{};
 app.renderMobilePanel=()=>{};
 app.handleAction('delete-leg',{side:'departure',id:selected});
 assert.equal(app.project.departureLegs.length,0);
 assert.equal(app.selected,null);
 app.handleAction('undo',{});
 assert.equal(app.project.departureLegs[0].id,selected);
});
