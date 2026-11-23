import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';

function fakeMap(){
 const objects=[];
 const create=(type,options)=>({type,options,handlers:{},on(name,fn){this.handlers[name]=fn;return this;},addTo(){objects.push(this);return this;},bindTooltip(){return this;}});
 const group=()=>({addTo(){return this;},clearLayers(){}});
 const map={setView(){return this;},on(){return this;},invalidateSize(){},getContainer(){return {style:{}};}};
 globalThis.L={map(){return map;},svg(){return {type:'renderer'};},tileLayer(){return {addTo(){}};},layerGroup:group,
    polyline(_pts,options){return create('polyline',options);},marker(_p,options){return create('marker',options);},circle(_p,options){return create('circle',options);},circleMarker(_p,options){return create('circleMarker',options);},divIcon(opts){return opts;}};
 return objects;
}
test('compact map exposes two endpoint plus controls and opens matching editors',async()=>{
 const oldWindow=globalThis.window;
 globalThis.window={matchMedia(){return {matches:true}}};
 const objects=fakeMap();const selected=[],points=[],added=[];
 try{
  const {PlannerMap}=await import('../../assets/js/ui/map.js');
  const p=createDemo(),planner=new PlannerMap({innerHTML:''},{onPick(){},onLegLength(){},onFinalLength(){},onHandle(){},onSelectLeg:(id)=>selected.push(id),onSelectPoint:(kind)=>points.push(kind),onAddAtEnd:(side)=>added.push(side)});
  planner.redraw(p,buildPlan(p));
  const plus=objects.filter(x=>x.type==='marker'&&x.options.icon.className==='glp-map-add-point');
  assert.equal(plus.length,2);
  plus[0].handlers.click();plus[1].handlers.click();assert.deepEqual(added,['departure','arrival']);
  const departureMarker=objects.find(x=>x.type==='marker'&&x.options.title.startsWith('離脱地点'));
  const landingMarker=objects.find(x=>x.type==='marker'&&x.options.title.startsWith('接地点'));
  departureMarker.handlers.click();landingMarker.handlers.click();assert.deepEqual(points,['departure','landing']);
  const hitLine=objects.find(x=>x.type==='polyline'&&x.options.weight===22);
  assert.ok(hitLine,'wide hit target must exist in compact mode');
  hitLine.handlers.click({latlng:{lat:35.4,lng:136.7}});
  assert.equal(selected.at(-1),p.departureLegs[0].id);
 }finally{globalThis.window=oldWindow;delete globalThis.L;}
});
test('desktop map does not overlay mobile add controls and hit areas',async()=>{
 const oldWindow=globalThis.window;globalThis.window={matchMedia(){return {matches:false}}};
 const objects=fakeMap();
 try{
  const {PlannerMap}=await import('../../assets/js/ui/map.js');
  const p=createDemo(),planner=new PlannerMap({innerHTML:''},{onPick(){},onLegLength(){},onFinalLength(){},onHandle(){},onAddAtEnd(){},onSelectLeg(){}});
  planner.redraw(p,buildPlan(p));
  assert.equal(objects.filter(x=>x.type==='marker'&&x.options.icon.className==='glp-map-add-point').length,0);
  assert.equal(objects.filter(x=>x.type==='polyline'&&x.options.weight===22).length,0);
 }finally{globalThis.window=oldWindow;delete globalThis.L;}
});
