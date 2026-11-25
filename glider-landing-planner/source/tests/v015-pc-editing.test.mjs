import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo, makeTurn, ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {parseProject, serializeProject} from '../../assets/js/io/project.js';

const priorDocument = globalThis.document;
globalThis.document = {querySelectorAll: () => []};
const {App} = await import('../../assets/js/main.js');
const {PlannerMap} = await import('../../assets/js/ui/map.js');
globalThis.document = priorDocument;

function fakeLeaflet() {
  const layers=[];
  const popupCalls=[];
  const map={
    on(){return this;},setView(){return this;},invalidateSize(){},
    closePopup(){popupCalls.push('close')},
    getContainer(){return {style:{}}}
  };
  const layer=(kind, coords, options={})=>({kind,coords,options,handlers:{},on(name,fn){this.handlers[name]=fn;return this;},
    addTo(){layers.push(this);return this;},bindTooltip(){return this;}});
  globalThis.L={
    map(){return map},svg(){return {}},tileLayer(){return {addTo(){}}},
    layerGroup(){return {addTo(){return this},clearLayers(){}}},
    polyline(pts,opts){return layer('polyline',pts,opts)},
    marker(pts,opts){return layer('marker',pts,opts)},
    circle(pts,opts){return layer('circle',pts,opts)},
    circleMarker(pts,opts){return layer('circleMarker',pts,opts)},
    divIcon(opts){return opts},
    popup(opts){
      const popup={options:opts,setLatLng(value){this.latlng=value;return this},setContent(value){this.content=value;return this},openOn(){popupCalls.push(this);return this}};
      return popup;
    }
  };
  return {layers,popupCalls,map};
}
function fakeElement(tag) {
  return {tag,children:[],handlers:{},appendChild(child){this.children.push(child)},addEventListener(name,fn){this.handlers[name]=fn}};
}

test('v0.1.15 desktop map: two smaller plus markers open menus and add the requested leg type', () => {
  const originalWindow=globalThis.window,originalDoc=globalThis.document,originalL=globalThis.L;
  globalThis.window={matchMedia:()=>({matches:false})};
  globalThis.document={createElement:fakeElement};
  try {
    const {layers,popupCalls}=fakeLeaflet(),calls=[];
    const project=previouslyPlannedProject();project.arrivalLegs=[makeTurn()];
    const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onAddAtEnd:(side,type)=>calls.push([side,type])});
    map.redraw(project,buildPlan(project));
    const plus=layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-add-point');
    assert.equal(plus.length,2);
    assert.deepEqual(plus.map(l=>l.options.icon.iconSize),[[36,36],[36,36]]);
    plus[0].handlers.click();
    const popup=popupCalls.at(-1);
    assert.match(popup.content.innerHTML,/離脱側へレグを追加/);
    assert.deepEqual(popup.content.children.map(b=>b.textContent),['＋ 直線','＋ 旋回']);
    popup.content.children[1].handlers.click();
    assert.equal(popupCalls.at(-1),'close');
    assert.deepEqual(calls,[['departure','turn']]);
    plus[1].handlers.click();
    popupCalls.at(-1).content.children[0].handlers.click();
    assert.deepEqual(calls,[['departure','turn'],['arrival','straight']]);
  } finally {globalThis.window=originalWindow;globalThis.document=originalDoc;globalThis.L=originalL;}
});

test('v0.1.15 compact map continues to use the bottom-sheet plus menu', () => {
  const originalWindow=globalThis.window,originalL=globalThis.L;
  globalThis.window={matchMedia:()=>({matches:true})};
  try {
    const {layers}=fakeLeaflet(),calls=[];
    const p=previouslyPlannedProject();p.arrivalLegs=[makeTurn()];const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onAddAtEnd:(side,type)=>calls.push([side,type])});
    map.redraw(p,buildPlan(p));
    const plus=layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-add-point');
    assert.equal(plus.length,2);
    assert.deepEqual(plus.map(l=>l.options.icon.iconSize),[[44,44],[44,44]]);
    plus[1].handlers.click();
    assert.deepEqual(calls,[['arrival',undefined]]);
  } finally {globalThis.window=originalWindow;globalThis.L=originalL;}
});

test('v0.1.15 point/length/final drag handles are absent while Bezier handle drag remains',()=>{
  const originalWindow=globalThis.window,originalL=globalThis.L;
  globalThis.window={matchMedia:()=>({matches:false})};
  try {
    const {layers}=fakeLeaflet(),selected=[];
    const p=previouslyPlannedProject();p.display.showHandles=true;
    const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onSelectPoint:x=>selected.push(x)});
    map.redraw(p,buildPlan(p));
    const markers=layers.filter(l=>l.kind==='marker');
    const pointMarkers=markers.filter(l=>l.options.title?.startsWith('離脱地点') || l.options.title?.startsWith('着陸地点'));
    assert.equal(pointMarkers.length,2);
    assert.ok(pointMarkers.every(m=>m.options.draggable !== true && !m.handlers.dragend));
    pointMarkers[0].handlers.click();pointMarkers[1].handlers.click();
    assert.deepEqual(selected,['departure','landing']);
    assert.equal(markers.filter(m=>m.options.icon.className==='glp-map-point').length,4,'two point markers plus two Bezier handles');
    assert.ok(markers.every(m=>!m.options.title?.includes('ドラッグで直線')&&!m.options.title?.includes('ファイナル長を変更')));
    const handles=markers.filter(m=>m.options.title==='Bezier制御ハンドル');
    assert.equal(handles.length,2);
    assert.ok(handles.every(m=>m.options.draggable===true && typeof m.handlers.dragend==='function'));
  } finally {globalThis.window=originalWindow;globalThis.L=originalL;}
});

test('v0.1.15 desktop sidebar scroll aligns offscreen editor without scrolling WordPress page',()=>{
  const originalWindow=globalThis.window;
  let pageScrolls=0;
  globalThis.window={matchMedia:()=>({matches:false}),scrollTo(){pageScrolls++}};
  try {
    const sidebar={scrollTop:100,getBoundingClientRect(){return {top:40,bottom:400}},scrollTo(opts){this.lastOptions=opts;this.scrollTop=opts.top}};
    const panel={getBoundingClientRect(){return {top:540,bottom:800}}};
    const app=Object.create(App.prototype);
    app.root={querySelector(selector){return selector==='.glp-sidebar-scroll'?sidebar:null}};
    app.scrollSidebarTo(panel);
    assert.equal(sidebar.scrollTop,588,'100 + 540 - 40 - 12');
    assert.equal(sidebar.lastOptions.behavior,'smooth');
    assert.equal(pageScrolls,0);
    const sameScroll=sidebar.scrollTop;
    app.scrollSidebarTo({getBoundingClientRect(){return {top:60,bottom:220}}});
    assert.equal(sidebar.scrollTop,sameScroll,'already visible sections do not move');
  }finally{globalThis.window=originalWindow;}
});

test('v0.1.15 selected desktop leg has an explicit delete button; numeric edits do not shift sidebar',()=>{
  const originalWindow=globalThis.window;
  globalThis.window={matchMedia:()=>({matches:false})};
  try {
    const app=Object.create(App.prototype),p=previouslyPlannedProject();
    p.departureLegs.push(makeTurn());
    app.project=p;app.plan=buildPlan(p);app.selected=p.departureLegs[1].id;
    const html=app.legList('departure');
    assert.match(html,/glp-desktop-leg-actions/);
    assert.match(html,/<button class="glp-btn glp-btn-sub glp-desktop-remove" data-action="delete-leg"/);
    assert.match(html,/このレグを削除/);
    const sidebar={scrollTop:235,innerHTML:''};
    app.root={querySelector(sel){return sel==='.glp-sidebar-scroll'?sidebar:null}};
    app.renderSidebar();
    assert.equal(sidebar.scrollTop,235,'normal re-render preserves scroll position');
  }finally{globalThis.window=originalWindow;}
});

test('v0.1.15 marker/route selection scrolls only the desktop editor',()=>{
  const originalWindow=globalThis.window;
  globalThis.window={matchMedia:()=>({matches:false})};
  try {
    const p=previouslyPlannedProject();
    const panel={dataset:{legId:p.departureLegs[0].id},getBoundingClientRect(){return {top:460,bottom:650}}};
    const sidebar={scrollTop:0,getBoundingClientRect(){return {top:40,bottom:360}},scrollTo({top}){this.scrollTop=top}};
    const app=Object.create(App.prototype);
    app.project=p;app.plan=buildPlan(p);app.selected=null;
    app.renderSidebar=()=>{};app.map={redraw(){}};app.updateCheckpointSelection=()=>{};
    app.root={querySelector(sel){return sel==='.glp-sidebar-scroll'?sidebar:null},querySelectorAll(sel){return sel==='.glp-leg[data-leg-id]'?[panel]:[]}};
    app.selectLeg(p.departureLegs[0].id,true);
    assert.equal(app.selected,p.departureLegs[0].id);
    assert.equal(sidebar.scrollTop,408);
  }finally{globalThis.window=originalWindow;}
});

test('current schema round-trips and map/print editing restrictions remain',()=>{
  const p=createDemo();
  assert.equal(ENGINE_VERSION,'0.1.17-test');
  assert.equal(p.schemaVersion,'1.3');
  assert.deepEqual(parseProject(serializeProject(p)),p);
  assert.deepEqual(p.performance,{bestLD:34,bestLDAtMps:25});
  const src=readFileSync(new URL('../src/ui/map.ts',import.meta.url),'utf8');
  assert.doesNotMatch(src,/onLegLength|onFinalLength|title: 'ドラッグで直線/);
  assert.match(src,/onHandle\(end,/);
  const print=readFileSync(new URL('../../assets/css/print.css',import.meta.url),'utf8');
  assert.match(print,/\.leaflet-popup-pane\{display:none!important\}/);
});
