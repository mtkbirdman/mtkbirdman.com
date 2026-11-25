import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {requiresLegacyAApproval} from '../../assets/js/io/project.js';
import {saveLegacyBackup,latestLegacyBackup} from '../../assets/js/io/autosave.js';

globalThis.document={querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');

function mockedApp(){
  const app=Object.create(App.prototype);
  app.project=createDemo();app.plan=buildPlan(app.project);
  app.history=[];app.future=[];app.editRevision=0;app.pendingLegacy=null;
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

test('legacy model A held for approval cannot be modified or autosaved',()=>{
  const {app}=mockedApp();
  const legacy=createDemo();legacy.schemaVersion='1.0';
  legacy.performance={model:'fixed-ld',fixedLD:30,bestLD:34,bestLDAtMps:25};
  assert.equal(requiresLegacyAApproval(legacy),true);
  app.pendingLegacy=structuredClone(legacy);
  const current=structuredClone(app.project);
  assert.equal(app.edit(p=>{p.name='should be blocked'}),false);
  assert.deepEqual(app.project,current);
  assert.equal(app.notice.includes('旧モデルA'),true);
  const oldWindow=globalThis.window;
  try{
    globalThis.window={confirm:()=>false};
    app.handleAction('legacy-approve',{});
    assert.deepEqual(app.project,current);
    assert.deepEqual(app.pendingLegacy,legacy);
    globalThis.window={confirm:()=>true};
    app.handleAction('legacy-approve',{});
    assert.equal(app.pendingLegacy,null);
    assert.equal(app.project.schemaVersion,'1.2');
    assert.equal(app.project.performance.model,'quadratic-drag');
    assert.equal(app.project.performance.quadratic.bestLD,30);
    assert.equal(app.history.length,0);
  } finally {globalThis.window=oldWindow;}
});

function fakeIndexedDb(){
  const stores=new Map();
  const db={
    objectStoreNames:{contains:name=>stores.has(name)},
    createObjectStore(name){stores.set(name,new Map());},
    close(){},
    transaction(name){
      const tx={oncomplete:null,onerror:null,objectStore(){return {
        put(value,key){stores.get(name).set(key,structuredClone(value));queueMicrotask(()=>tx.oncomplete?.());},
        get(key){const request={};queueMicrotask(()=>{request.result=structuredClone(stores.get(name).get(key));request.onsuccess?.();});return request;}
      };}};
      return tx;
    }
  };
  return {open(){
    const req={result:db};
    queueMicrotask(()=>{req.onupgradeneeded?.();req.onsuccess?.();});
    return req;
  }};
}

test('IndexedDB legacy backup is independent from the rolling snapshot store and preserves exact JSON',async()=>{
  const oldWindow=globalThis.window,oldIndexedDb=globalThis.indexedDB;
  const old=createDemo();old.schemaVersion='1.0';old.performance={model:'fixed-ld',fixedLD:29,bestLD:34,bestLDAtMps:25};
  try{
    globalThis.window={indexedDB:{}};
    globalThis.indexedDB=fakeIndexedDb();
    await saveLegacyBackup(old);
    const backup=await latestLegacyBackup();
    assert.deepEqual(backup,old);
    assert.equal(backup.performance.model,'fixed-ld');
  } finally{globalThis.window=oldWindow;globalThis.indexedDB=oldIndexedDb;}
});
