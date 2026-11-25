import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';

globalThis.document={querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');

function mockedApp(){
  const app=Object.create(App.prototype);
  app.project=createDemo();app.plan=buildPlan(app.project);
  app.history=[];app.future=[];app.editRevision=0;
  app.polarEditorOpen=false;app.polarDraft='';app.selected=null;
  app.map={fit(){}};
  let saves=0;
  app.render=()=>{};app.renderResults=()=>{};app.renderSidebar=()=>{};
  app.scheduleSave=()=>{saves++};
  return {app,get saves(){return saves}};
}

test('C draft remains separate until Apply and a malformed draft does not overwrite valid B',()=>{
  const {app}=mockedApp();const before=structuredClone(app.project);
  app.polarEditorOpen=true;app.polarDraft='70,0.7\n80,0.6';
  app.handleAction('polar-apply',{});
  assert.deepEqual(app.project,before);
  assert.equal(app.polarEditorOpen,true);
  assert.match(app.notice,/4〜64/);
  app.polarDraft='100, 1.0\n70, 1.2\n90, 0.7\n80, 0.85';
  app.handleAction('polar-apply',{});
  assert.equal(app.project.performance.model,'polar-table');
  assert.equal(app.project.performance.polarTable.points.length,4);
  assert.ok(app.project.performance.polarTable.points.every((p,i,a)=>i===0||a[i-1].speedMps<p.speedMps));
  assert.equal(app.history.length,1);
  assert.equal(app.polarEditorOpen,false);
});

test('C draft Cancel leaves model and data untouched; B/C editing keeps old C points',()=>{
  const {app}=mockedApp();
  app.polarDraft='70,1\n80,.8\n90,.75\n100,1';
  app.handleAction('polar-apply',{});
  const before=structuredClone(app.project.performance);
  app.handleAction('polar-edit',{});
  assert.ok(app.polarDraft.includes('70'));
  app.polarDraft='not valid';app.handleAction('polar-cancel',{});
  assert.deepEqual(app.project.performance,before);
  app.edit(p=>{p.performance.model='quadratic-drag';});
  assert.deepEqual(app.project.performance.polarTable,before.polarTable);
  app.edit(p=>{p.performance.model='polar-table';});
  assert.deepEqual(app.project.performance,before);
});

test('valid B/C editing remains available without a legacy approval gate',()=>{
  const {app}=mockedApp();
  assert.equal(app.edit(p=>{p.name='編集済み計画'}),true);
  assert.equal(app.project.name,'編集済み計画');
  assert.equal(app.history.length,1);
});
