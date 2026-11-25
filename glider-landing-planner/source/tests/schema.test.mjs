import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, existsSync} from 'node:fs';
import {createDemo, makeStraight, ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {loadProject, parseProject, serializeProject, validateProject} from '../../assets/js/io/project.js';
import {reportHtml} from '../../assets/js/ui/print.js';

globalThis.document={querySelectorAll:()=>[],querySelector:()=>null};
const {App}=await import('../../assets/js/main.js');

test('all pre-1.3 schemas and unknown schemas are rejected without mutation',()=>{
  for(const version of ['1.0','1.1','1.2','9.9',undefined]){
    const p=createDemo();p.schemaVersion=version;
    const before=structuredClone(p);
    assert.throws(()=>loadProject(p),/schemaVersion/);
    assert.throws(()=>parseProject(JSON.stringify(p)),/schemaVersion/);
    assert.throws(()=>serializeProject(p),/schemaVersion/);
    assert.deepEqual(p,before);
  }
});

test('1.3 requires exactly the two valid performance parameters and does not fill absent fields',()=>{
  for(const performance of [{bestLD:34}, {bestLDAtMps:25}, {bestLD:1,bestLDAtMps:25}, {bestLD:34,bestLDAtMps:0}, {bestLD:NaN,bestLDAtMps:25}, {bestLD:34,bestLDAtMps:25,model:'quadratic-drag'}, {bestLD:34,bestLDAtMps:25,polarTable:null}]){
    const p=createDemo();p.performance=performance;
    assert.throws(()=>validateProject(p),/性能パラメータ/);
    assert.throws(()=>loadProject(p),/性能パラメータ/);
  }
  const p=createDemo();delete p.display.referenceLineOpacity;
  assert.throws(()=>loadProject(p),/不透明度/);
});

test('performance controls remain editable; no selector or table UI at either viewport size',()=>{
  const beforeWindow=globalThis.window, beforeInput=globalThis.HTMLInputElement;
  try {
    globalThis.HTMLInputElement=class {};
    for(const compact of [false,true]){
      globalThis.window={matchMedia:()=>({matches:compact})};
      const p=createDemo();p.departureLegs=[makeStraight(300)];
      const app=Object.create(App.prototype), sidebar={innerHTML:'',scrollTop:0};
      Object.assign(app,{project:p,plan:buildPlan(p),history:[],future:[],editRevision:0,selected:null,mobilePanel:null,root:{querySelector:()=>sidebar},render(){},renderResults(){},scheduleSave(){}});
      app.renderSidebar();
      assert.match(sidebar.innerHTML,/data-path="performance\.bestLD"/);
      assert.match(sidebar.innerHTML,/data-path="performance\.bestLDAtMps"/);
      assert.doesNotMatch(sidebar.innerHTML,/沈下計算モデル|ポーラー|performance\.model/);
      const initial=app.plan.departure.end.h;
      const input=Object.assign(new HTMLInputElement(),{dataset:{path:'performance.bestLD'},type:'number',value:'25'});
      app.handleChange(input);
      assert.equal(app.project.performance.bestLD,25);
      assert.ok(app.plan.departure.end.h<initial);
      input.dataset={path:'performance.bestLDAtMps',scale:'3.6'};input.value='108';
      app.handleChange(input);
      assert.equal(app.project.performance.bestLDAtMps,30);
      assert.deepEqual(parseProject(serializeProject(app.project)),app.project);
    }
  }finally{globalThis.window=beforeWindow;globalThis.HTMLInputElement=beforeInput;}
});

test('print HTML includes both performance values with no model selector labels',()=>{
  const p=createDemo();p.performance={bestLD:38,bestLDAtMps:30};
  const mapCopy={querySelectorAll:()=>[],removeAttribute(){},classList:{add(){}},style:{},outerHTML:'<div class="leaflet"></div>'};
  const html=reportHtml({project:p,plan:buildPlan(p),mapEl:{clientWidth:900,clientHeight:420,cloneNode:()=>mapCopy},metricsHtml:'',diagnosticsHtml:'',title:'GLP',engineVersion:ENGINE_VERSION});
  assert.match(html,/最大滑空比/);assert.match(html,/>38</);assert.match(html,/最良滑空速度/);assert.match(html,/108 km\/h/);
  assert.doesNotMatch(html,/沈下モデル|速度別ポーラー|B：簡易抗力/);
});

test('removed source and compiled modules and model-specific code are absent',()=>{
  for(const file of ['../src/flight/polar.ts','../../assets/js/flight/polar.js'])assert.equal(existsSync(new URL(file,import.meta.url)),false);
  for(const file of ['../src/main.ts','../src/flight/core.ts','../src/io/project.ts','../src/planning/connector.ts','../src/planning/plan.ts','../src/model/types.ts','../src/ui/print.ts','../../assets/css/app.css']){
    assert.doesNotMatch(readFileSync(new URL(file,import.meta.url),'utf8'),/polarTable|polar-table|PolarEvaluationError|POLAR_|pchip|performance\.quadratic|performance\.model|glp-polar/);
  }
});

test('README example works through the real edit handler and current JSON round-trip',()=>{
  const p=createDemo(),app=Object.create(App.prototype);
  Object.assign(app,{project:p,plan:buildPlan(p),history:[],future:[],editRevision:0,selected:null,mobilePanel:null,render(){},renderResults(){},scheduleSave(){}});
  const priorWindow=globalThis.window;globalThis.window={matchMedia:()=>({matches:false}),confirm:()=>true};
  try{
    app.handleAction('preset-left',{});
    assert.equal(app.project.arrivalLegs.length,4);
    for(let i=0;i<2;i++)assert.equal(app.edit(p=>p.departureLegs.push(makeStraight(300))),true);
    assert.equal(app.edit(p=>{p.defaultSpeedMps=95/3.6;}),true);
    assert.equal(app.plan.connector.feasible,true);
    assert.equal(app.plan.notices.filter(n=>n.level==='error').length,0);
    assert.deepEqual(parseProject(serializeProject(app.project)),app.project);
  }finally{globalThis.window=priorWindow;}
});
