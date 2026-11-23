import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo,trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {checkpointRows,checkpointTable,landingRelativeHeight} from '../../assets/js/ui/checkpoints.js';
import {heightForDisplay,heightFromDisplay,heightInputText} from '../../assets/js/ui/height-units.js';
import {serializeProject,parseProject} from '../../assets/js/io/project.js';
const near=(a,b,tol=1e-9)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);

test('v0.1.4 height display defaults to ft and converting input does not mutate canonical MSL meters',()=>{
 const p=createDemo();assert.equal(p.display.altitudeUnit,'ft');
 p.departure.altitudeMSLM=500;p.landing.elevationMSLM=130;
 assert.equal(heightInputText(p.departure.altitudeMSLM,'ft'),'1640.42');
 near(heightFromDisplay(heightForDisplay(500,'ft'),'ft'),500);
 near(heightFromDisplay(heightForDisplay(130,'m'),'m'),130);
 const original=serializeProject(p);p.display.altitudeUnit='m';
 assert.equal(heightInputText(p.departure.altitudeMSLM,'m'),'500');
 p.display.altitudeUnit='ft';
 assert.equal(heightInputText(p.landing.elevationMSLM,'ft'),'426.51');
 assert.equal(p.departure.altitudeMSLM,500);assert.equal(p.landing.elevationMSLM,130);
 const restored=parseProject(original);assert.equal(restored.departure.altitudeMSLM,500);
});

test('v0.1.4 checkpoint table uses landing-relative height, shows 0 at landing and no calculation-type column',()=>{
 const p=createDemo();p.landing.elevationMSLM=130;p.departure.altitudeMSLM=500;
 p.arrivalLegs=trafficPreset('right');
 const plan=buildPlan(p),rows=checkpointRows(plan,p);
 const arrivalTarget=plan.checkpoints.find(c=>c.kind==='landing');assert.ok(arrivalTarget);
 near(landingRelativeHeight(arrivalTarget,p),0);
 const departure=plan.checkpoints[0];near(landingRelativeHeight(departure,p),370);
 const i=rows.findIndex(r=>r.label==='ダウンウィンド開始');assert.ok(i>0);
 assert.equal(rows[i-1].kind,'margin');
 assert.equal(rows[i-1].label,'ダウンウィンド開始の高度余裕');
 const table=checkpointTable(plan,p);
 assert.match(table,/滑空場基準高度（ft）/);
 assert.match(table,/MSL高度（ft）/);
 assert.doesNotMatch(table,/計算の種類/);
 assert.doesNotMatch(table,/<tr class="glp-cp-group"/);
 assert.doesNotMatch(table,/離脱側（予想高度）/);
 assert.doesNotMatch(table,/着陸側（逆算必要高度）/);
 assert.match(table,/進入目標点<\/td><td>0 ft<\/td><td>427 ft<\/td>/);
 assert.match(table,/ダウンウィンド開始の高度余裕<\/td><td>[^<]+<\/td><td>—<\/td>/);
 p.display.altitudeUnit='m';
 assert.match(checkpointTable(plan,p),/滑空場基準高度（m）/);
 assert.match(checkpointTable(plan,p),/離脱地点<\/td><td>370 m<\/td><td>500 m<\/td>/);
});

test('v0.1.4 legacy saved JSON still imports and retains original MSL fields',()=>{
 const p=createDemo();p.engineVersion='0.1.2-test';p.departure.altitudeMSLM=520.1;
 p.landing.elevationMSLM=35.5;p.display.altitudeUnit='m';
 const imported=parseProject(serializeProject(p));
 assert.equal(imported.engineVersion,'0.1.8-test');
 near(imported.departure.altitudeMSLM,520.1);near(imported.landing.elevationMSLM,35.5);
 assert.equal(imported.display.altitudeUnit,'m');
});
