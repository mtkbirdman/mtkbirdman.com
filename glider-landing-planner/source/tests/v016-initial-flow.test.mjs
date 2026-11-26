import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo, makeStraight, makeTurn, trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {LocalFrame,headingVector,dot,sub} from '../../assets/js/geometry/geo.js';
import {chartSvg} from '../../assets/js/ui/chart.js';
import {parseProject,serializeProject} from '../../assets/js/io/project.js';
import {previouslyPlannedProject} from './previously-planned-fixture.mjs';

const existingDocument=globalThis.document;
globalThis.document={querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');
const {PlannerMap}=await import('../../assets/js/ui/map.js');
globalThis.document=existingDocument;

function leaflets(){
 const layers=[],popups=[];
 const layer=(kind,coords,options={})=>({kind,coords,options,handlers:{},on(name,cb){this.handlers[name]=cb;return this;},addTo(){layers.push(this);return this;},bindTooltip(){return this;}});
 const map={on(){return this;},setView(){return this;},invalidateSize(){},closePopup(){popups.push('close');},getContainer(){return {style:{}}}};
 globalThis.L={map(){return map},svg(){return {}},tileLayer(){return {addTo(){return this}}},layerGroup(){return {addTo(){return this},clearLayers(){}}},
   polyline(coords,opts){return layer('line',coords,opts)},circle(coords,opts){return layer('circle',coords,opts)},circleMarker(coords,opts){return layer('point',coords,opts)},marker(coords,opts){return layer('marker',coords,opts)},divIcon(opts){return opts},
   popup(options){return {options,setLatLng(value){this.at=value;return this},setContent(value){this.content=value;return this},openOn(){popups.push(this);return this}}}};
 return {layers,popups};
}
function fakeElement(tag){return {tag,innerHTML:'',textContent:'',children:[],events:{},appendChild(c){this.children.push(c)},addEventListener(event,cb){this.events[event]=cb}};}
function withApp(project=createDemo()){
 const app=Object.create(App.prototype);
 app.project=project;app.plan=buildPlan(project);
 app.history=[];app.future=[];app.editRevision=0;app.notice='';app.status='';app.selected=null;app.mobilePanel=null;
 app.render=()=>{};app.scheduleSave=()=>{};
 return app;
}

test('new plan is unplanned, with no synthetic flight leg, Bezier, or false connection error',()=>{
 const p=createDemo(),r=buildPlan(p);
 assert.equal(p.schemaVersion,'1.3');assert.equal(p.departureLegs.length,0);assert.equal(p.arrivalLegs.length,0);
 assert.equal(p.connector,null);assert.equal(r.connector,null);assert.equal(r.departure.segments.length,0);
 assert.ok(r.notices.some(n=>n.code==='PLAN_NOT_STARTED'&&n.level==='info'));
 assert.ok(!r.notices.some(n=>n.code==='CONNECTOR_MISSING'));
 assert.equal(r.altitudeValid,false);
 assert.equal(r.forwardAll.length,1);
 assert.match(chartSvg(r,p),/着陸側レグまたは場周プリセットの追加後/);
 assert.doesNotMatch(chartSvg(r,p),/glp-departure-altitude/);
});

test('departure directional guide is exactly 250 m from the anchor and does not enter the plan',()=>{
 const previousL=globalThis.L,previousW=globalThis.window;
 globalThis.window={matchMedia:()=>({matches:false})};
 try{
   const {layers}=leaflets(),p=createDemo();
   const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onSelectPoint(){},onAddAtEnd(){}});
   map.redraw(p,buildPlan(p));
   const guides=layers.filter(l=>l.kind==='line' && l.options.color==='#2563eb' && l.options.dashArray==='7,6');
   assert.equal(guides.length,1);
   assert.equal(guides[0].options.interactive,false,'heading guide cannot intercept map clicks');
   const frame=new LocalFrame(p.landing.point);
   const from=frame.toLocal({latDeg:guides[0].coords[0][0],lonDeg:guides[0].coords[0][1]});
   const to=frame.toLocal({latDeg:guides[0].coords[1][0],lonDeg:guides[0].coords[1][1]});
   const displacement=sub(to,from),heading=headingVector(95);
   assert.ok(Math.abs(Math.hypot(displacement.e,displacement.n)-250)<.001);
   assert.ok(Math.abs(dot(displacement,heading)-250)<.001);
   assert.equal(layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-add-point').length,0);
   assert.equal(layers.filter(l=>l.kind==='point').length,0,'unplanned checkpoints do not clutter the map');
   assert.ok(layers.some(l=>l.kind==='line'&&l.options.color==='#059669'&&l.options.dashArray==='9,4'),'landing direction guide');
   assert.ok(layers.some(l=>l.kind==='line'&&l.options.color==='#059669'&&l.options.weight===2),'landing arrowhead');
 }finally{globalThis.L=previousL;globalThis.window=previousW;}
});

test('desktop initial departure marker offers add-straight, add-turn, and edit; landing opens preset editor',()=>{
 const prevL=globalThis.L,prevW=globalThis.window,prevD=globalThis.document;
 globalThis.window={matchMedia:()=>({matches:false})};globalThis.document={createElement:fakeElement};
 try{
   const {layers,popups}=leaflets(),added=[],selected=[];
   const p=createDemo(),map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onSelectPoint:x=>selected.push(x),onAddAtEnd:(s,t)=>added.push([s,t])});
   map.redraw(p,buildPlan(p));
   const pointMarkers=layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-point');
   assert.equal(pointMarkers.length,2);
   pointMarkers[0].handlers.click();
   const popup=popups.at(-1);
   assert.deepEqual(popup.content.children.map(c=>c.textContent),['＋ 直線','＋ 旋回','離脱地点・条件を編集']);
   popup.content.children[1].events.click();
   assert.deepEqual(added,[['departure','turn']]);
   pointMarkers[0].handlers.click();popups.at(-1).content.children[2].events.click();
   pointMarkers[1].handlers.click();
   assert.deepEqual(selected,['departure','landing']);
 }finally{globalThis.L=prevL;globalThis.window=prevW;globalThis.document=prevD;}
});

test('compact initial departure marker opens addition panel; landing marker opens preset/landing editor',()=>{
 const prevL=globalThis.L,prevW=globalThis.window;
 globalThis.window={matchMedia:()=>({matches:true})};
 try{
   const {layers}=leaflets(),added=[],selected=[];
   const p=createDemo(),map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onSelectPoint:x=>selected.push(x),onAddAtEnd:(s,t)=>added.push([s,t])});
   map.redraw(p,buildPlan(p));
   const pointMarkers=layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-point');
   pointMarkers[0].handlers.click();pointMarkers[1].handlers.click();
   assert.deepEqual(added,[['departure',undefined]]);
   assert.deepEqual(selected,['landing']);
   const app=withApp();app.mobilePanel={kind:'add',side:'departure'};
   let content='';const panel={hidden:true,querySelector(sel){if(sel==='.glp-mobile-editor-head')return {innerHTML:''};if(sel==='.glp-mobile-editor-content')return {get scrollTop(){return 0},set scrollTop(_){},set innerHTML(value){content=value}};}};
   app.root={querySelector(sel){return sel==='.glp-mobile-editor'?panel:null}};
   app.renderMobilePanel();
   assert.match(content,/data-action="add-straight"/);assert.match(content,/data-action="add-turn"/);assert.match(content,/data-action="edit-departure"/);
 }finally{globalThis.L=prevL;globalThis.window=prevW;}
});

test('first departure leg removes guide and adds only departure extension button',()=>{
 const prevL=globalThis.L,prevW=globalThis.window;
 globalThis.window={matchMedia:()=>({matches:false})};
 try{
   const {layers}=leaflets(),p=createDemo();p.departureLegs=[makeStraight(300)];
   const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){},onAddAtEnd(){}});
   map.redraw(p,buildPlan(p));
   assert.equal(layers.filter(l=>l.kind==='line'&&l.options.dashArray==='7,6'&&l.options.color==='#2563eb').length,0);
   const plus=layers.filter(l=>l.kind==='marker'&&l.options.icon.className==='glp-map-add-point');
   assert.equal(plus.length,1);assert.match(plus[0].options.title,/離脱側/);
 }finally{globalThis.L=prevL;globalThis.window=prevW;}
});

test('arrival first leg automatically creates and optimizes Bezier; deleting last arrival leg restores unplanned state',()=>{
 const app=withApp();
 assert.equal(app.edit(p=>p.departureLegs.push(makeTurn())),true);
 assert.equal(app.project.connector,null);
 assert.equal(app.edit(p=>p.arrivalLegs.unshift(makeStraight(300))),true);
 assert.ok(app.project.connector);
 assert.equal(app.project.connector.adjustmentMode,'auto');
 assert.ok(app.plan.connector?.points.length>2);
 const saved=structuredClone(app.project);
 assert.equal(app.edit(p=>p.arrivalLegs.splice(0,1)),true);
 assert.equal(app.project.connector,null);
 assert.ok(app.plan.notices.some(n=>n.code==='PLAN_NOT_STARTED'));
 assert.equal(app.history.length,3);
 app.future.push(structuredClone(app.project));app.project=app.history.pop();app.plan=buildPlan(app.project);
 assert.deepEqual(app.project,saved,'undo restores both the arrival leg and its resolved connector');
});

test('preset action restores a connector after a manually disconnected populated arrival',()=>{
 const priorWindow=globalThis.window;
 globalThis.window={matchMedia:()=>({matches:false}),confirm:()=>true};
 try{
   const app=withApp();
   app.handleAction('preset-left',{});
   assert.ok(app.project.connector);
   assert.equal(app.project.arrivalLegs.length,4);
   app.edit(p=>p.connector=null);
   assert.equal(app.project.connector,null);
   app.handleAction('preset-right',{});
   assert.ok(app.project.connector);
   assert.equal(app.project.connector.adjustmentMode,'auto');
   assert.equal(app.project.arrivalLegs[0].checkpointLabel,'ダウンウィンド開始');
 }finally{globalThis.window=priorWindow;}
});

test('saved current-schema project retains its existing legs and connector, while new JSON round-trips unplanned data',()=>{
 const old=previouslyPlannedProject();old.arrivalLegs=trafficPreset('left');
 const original=structuredClone(old);
 const loaded=parseProject(serializeProject(old));
 assert.deepEqual(loaded,original);
 const p=createDemo();assert.deepEqual(parseProject(serializeProject(p)),p);
 const source=readFileSync(new URL('../src/main.ts',import.meta.url),'utf8');
 assert.match(source,/data-point-editor="landing"/);
 assert.match(source,/着陸地点・ファイナル・場周設定/);
 assert.match(source,/右場周プリセット/);
});

test('version and shortcode stay aligned for FF FTP release',()=>{
 const php=readFileSync(new URL('../../glider-landing-planner.php',import.meta.url),'utf8');
 assert.match(php,/Version: 0\.1\.20/);
 assert.match(php,/add_shortcode\('glider_landing_planner'/);
 const changelog=readFileSync(new URL('../../CHANGELOG.md',import.meta.url),'utf8');
 assert.match(changelog,/## \[0\.1\.20\]/);
});


test('unplanned result cards and diagnostic show no fabricated flight duration, distance, or margin',()=>{
 const app=withApp(),elements=new Map();
 app.root={querySelector(selector){if(!elements.has(selector))elements.set(selector,{innerHTML:''});return elements.get(selector);},querySelectorAll(){return []}};
 app.renderResults();
 const metrics=elements.get('.glp-metrics').innerHTML;
 assert.match(metrics,/推定総距離<\/small><strong>—<\/strong>/);
 assert.match(metrics,/推定飛行時間<\/small><strong>—<\/strong>/);
 assert.match(metrics,/高度余裕<\/small><strong>—<\/strong>/);
 assert.match(elements.get('.glp-diagnostics').innerHTML,/未設計/);
});

test('real preset action opens connector; deleting last arrival leg and Undo restore the state',()=>{
 const prevWindow=globalThis.window;
 globalThis.window={matchMedia:()=>({matches:false}),confirm:()=>true};
 try{
   const app=withApp();
   app.renderSidebar=()=>{};app.map={redraw(){}};app.updateCheckpointSelection=()=>{};
   app.root={querySelector(){return null},querySelectorAll(){return []}};
   app.handleAction('preset-left',{});
   assert.equal(app.project.arrivalLegs.length,4);
   assert.ok(app.project.connector);
   const saved=structuredClone(app.project);
   // The preset is four legs: delete from the end until the last one is gone.
   for(const leg of [...app.project.arrivalLegs]) app.handleAction('delete-leg',{side:'arrival',id:leg.id});
   assert.equal(app.project.arrivalLegs.length,0);
   assert.equal(app.project.connector,null);
   app.handleAction('undo',{});
   assert.equal(app.project.arrivalLegs.length,1);
   assert.ok(app.project.connector,'Undo restores the connector on the last leg');
   for(let i=0;i<3;i++) app.handleAction('undo',{});
   assert.deepEqual(app.project,saved);
 }finally{globalThis.window=prevWindow;}
});

test('manual addition on either side defaults to ground-track without rewriting saved legs',()=>{
 const p=createDemo();p.engineVersion='0.1.17-test';p.wind={fromDeg:270,speedMps:5};
 p.departureLegs=[{...makeStraight(),guidance:'heading'}];
 const old=parseProject(serializeProject(p));assert.deepEqual(old,p);
 const app=withApp(old);app.selectLeg=()=>{};
 for(const side of ['departure','arrival']) for(const type of ['straight','turn']){
   app.handleAction(`add-${type}`,{side});
   const added=side==='departure'?app.project.departureLegs.at(-1):app.project.arrivalLegs[0];
   assert.equal(added.type,type);assert.equal(added.guidance,'ground-track');
   assert.equal(added.directionDeg,null);
   if(type==='turn')assert.equal(added.angleInputMode,'angle');
 }
 assert.deepEqual(app.project.departureLegs[0],p.departureLegs[0]);
 assert.deepEqual(parseProject(serializeProject(app.project)),app.project);
});
