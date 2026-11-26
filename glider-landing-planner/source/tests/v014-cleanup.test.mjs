import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo, makeTurn} from '../../assets/js/model/defaults.js';
import * as projectIO from '../../assets/js/io/project.js';
import * as autosave from '../../assets/js/io/autosave.js';
import {buildPlan} from '../../assets/js/planning/plan.js';

// Old A is no longer a supported data format: no approval, backup or conversion API remains.
test('removed model A APIs are absent; current project still round-trips',()=>{
  assert.equal('requiresLegacyAApproval' in projectIO,false);
  assert.equal('saveLegacyBackup' in autosave,false);
  assert.equal('latestLegacyBackup' in autosave,false);
  const project=createDemo();
  assert.deepEqual(projectIO.parseProject(projectIO.serializeProject(project)),project);
});

test('unsupported schema 1.0 model cannot be auto-converted and input remains unchanged',()=>{
  const input=createDemo();
  input.schemaVersion='1.0';
  input.performance={model:'fixed-ld',fixedLD:30,bestLD:34,bestLDAtMps:25};
  const original=structuredClone(input);
  assert.throws(()=>projectIO.loadProject(input),/schemaVersion/);
  assert.deepEqual(input,original);
});


test('turn editor removes auto-calculation label but keeps readonly and target editing',async()=>{
  globalThis.document={querySelectorAll:()=>[]};
  const {App}=await import('../../assets/js/main.js');
  const app=Object.create(App.prototype);
  app.project=createDemo();
  app.project.departureLegs=[{...makeTurn('right',450),angleInputMode:'target',targetDirectionDeg:180}];
  app.plan=buildPlan(app.project);
  const targetHtml=app.legEditor('departure',app.project.departureLegs[0],0);
  assert.match(targetHtml, /進行方向変化量 <em>°<\/em><\/span><input type="number" readonly/);
  assert.match(targetHtml,/data-path="departureLegs\.0\.targetDirectionDeg"/);
  assert.doesNotMatch(targetHtml,/（自動計算）/);
  app.project.departureLegs[0].angleInputMode='angle';
  app.project.departureLegs[0].targetDirectionDeg=null;
  app.plan=buildPlan(app.project);
  const angleHtml=app.legEditor('departure',app.project.departureLegs[0],0);
  assert.match(angleHtml,/data-path="departureLegs\.0\.angleDeg"/);
  assert.match(angleHtml,/旋回後の進行方向 <em>°<\/em><\/span><input type="number" readonly/);
  assert.doesNotMatch(angleHtml,/（自動計算）/);
});

test('no model-A migration controls or autosave actions remain in application source',()=>{
  const main=readFileSync(new URL('../src/main.ts',import.meta.url),'utf8');
  const io=readFileSync(new URL('../src/io/autosave.ts',import.meta.url),'utf8');
  assert.doesNotMatch(main,/pendingLegacy|legacy-approve|legacy-backup|requiresLegacyAApproval|saveLegacyBackup/);
  assert.doesNotMatch(io,/legacy-backups|saveLegacyBackup|latestLegacyBackup/);
});
