import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo,makeTurn,trafficPreset,ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {parseProject,serializeProject,validateProject} from '../../assets/js/io/project.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {checkpointTable,checkpointRows} from '../../assets/js/ui/checkpoints.js';

const read=(relative)=>readFileSync(new URL(relative,import.meta.url),'utf8');

test('v0.1.4: new plan and newly added turns start at 30 degrees, handles default off',()=>{
  const p=createDemo();
  assert.equal(ENGINE_VERSION,'0.1.16-test');
  assert.equal(p.engineVersion,ENGINE_VERSION);
  assert.equal(p.name,'新規計画');
  assert.equal(p.display.showHandles,false);
  assert.equal(p.display.altitudeUnit,'ft');
  assert.equal(makeTurn().bankDeg,30);
  for(const s of ['left','right']) assert.ok(trafficPreset(s).filter(l=>l.type==='turn').every(l=>l.bankDeg===30));
});

test('v0.1.4: existing legacy departure/arrival turn banks migrate once to 30 and keep other fields',()=>{
  const p=createDemo();
  p.engineVersion='0.1.3-test';
  p.name='新規計画（架空の初期値）';
  p.departureLegs.push({...makeTurn('left',270),bankDeg:25});
  p.arrivalLegs=trafficPreset('right');
  p.arrivalLegs[1].bankDeg=18;
  p.arrivalLegs[3].bankDeg=42;
  p.display.showHandles=true;
  const migrated=parseProject(JSON.stringify(p));
  assert.equal(migrated.engineVersion,ENGINE_VERSION);
  assert.equal(migrated.name,'新規計画');
  assert.deepEqual([...migrated.departureLegs,...migrated.arrivalLegs].filter(l=>l.type==='turn').map(l=>l.bankDeg),[30,30,30]);
  assert.equal(migrated.departureLegs.at(-1).angleDeg,270);
  assert.equal(migrated.display.showHandles,true);
  // Subsequent edits in the NEW engine version must NOT be silently rewritten.
  migrated.departureLegs.at(-1).bankDeg=37;
  migrated.arrivalLegs[1].bankDeg=29;
  assert.equal(validateProject(migrated).departureLegs.at(-1).bankDeg,37);
  const again=parseProject(serializeProject(migrated));
  assert.equal(again.departureLegs.at(-1).bankDeg,37);
  assert.equal(again.arrivalLegs[1].bankDeg,29);
});

test('v0.1.4: auto-save recovery is silent and guarded against overwriting user edits',()=>{
  const main=read('../src/main.ts');
  assert.doesNotMatch(main,/window\.confirm\('前回の自動保存プロジェクトがあります/);
  assert.match(main,/this\.editRevision\s*!==\s*revisionAtStartup/);
  assert.match(main,/this\.status\s*=\s*'自動保存から復元しました'/);
  assert.match(main,/自動保存データが不正なため復元できません/);
  assert.match(main,/void saveSnapshot\(this\.project\)/);
  assert.match(main,/\[\['ft',\s*'ft'\],\s*\['m',\s*'m'\]\]/);
});

test('v0.1.4: height table shows independent airfield-relative and MSL values without group rows',()=>{
  const p=createDemo();
  p.departure.altitudeMSLM=500;
  p.landing.elevationMSLM=130;
  p.arrivalLegs=trafficPreset('left');
  const plan=buildPlan(p);
  const rows=checkpointRows(plan,p);
  assert.ok(rows.every(r=>r.kind==='height'||r.kind==='margin'));
  const marginIndex=rows.findIndex(r=>r.kind==='margin');
  assert.equal(rows[marginIndex+1].label,'ダウンウィンド開始');
  const landing=rows.find(r=>r.label==='進入目標点');
  assert.equal(landing.relativeM,0);
  assert.equal(landing.mslM,130);
  const table=checkpointTable(plan,p);
  assert.match(table,/滑空場基準高度（ft）/);
  assert.match(table,/MSL高度（ft）/);
  assert.match(table,/進入目標点<\/td><td>0 ft<\/td><td>427 ft/);
  assert.match(table,/ダウンウィンド開始の高度余裕<\/td><td>[^<]+<\/td><td>—/);
  assert.doesNotMatch(table,/glp-cp-group/);
  p.display.altitudeUnit='m';
  assert.match(checkpointTable(plan,p),/進入目標点<\/td><td>0 m<\/td><td>130 m/);
  assert.match(checkpointTable(plan,p,true),/class="glp-print-table"/);
});

test('v0.1.4: left menu no fixed height, table no vertical inner scroll',()=>{
  const css=read('../../assets/css/app.css');
  assert.match(css,/\.glp-sidebar-scroll\{position:sticky;top:40px;max-height:calc\(100dvh - 52px\)/);
  assert.match(css,/\.glp-tablewrap\{max-height:none;overflow-x:auto;overflow-y:visible/);
  assert.doesNotMatch(css,/height:min\(1050px/);
  assert.match(read('../../glider-landing-planner.php'),/Version: 0\.1\.16/);
});
