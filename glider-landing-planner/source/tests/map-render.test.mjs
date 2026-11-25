import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo,trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';

// Renderer regression test: overlays must be vector SVG, never theme-sensitive Canvas.
test('GLP WordPress/Cocoon vector paths use explicit SVG renderer', async () => {
  const recorded={mapOptions:null,renderers:[],paths:[],circles:[],markers:[]};
  const makeLayer=(kind,record,options={})=>({
    kind,options,
    addTo(group){record.push(this);return this;},
    on(){return this;},bindTooltip(){return this;},
  });
  const makeGroup=()=>({items:[],addTo(){return this;},clearLayers(){this.items=[];}});
  globalThis.L={
    map(_element,options){recorded.mapOptions=options;return {setView(){return this;},on(){return this;},invalidateSize(){},fitBounds(){}};},
    svg(options){const renderer={svgRenderer:true,options};recorded.renderers.push(renderer);return renderer;},
    tileLayer(){return {addTo(){return this;}};},
    layerGroup:makeGroup,
    polyline(_points,options){return makeLayer('polyline',recorded.paths,options);},
    circle(_points,options){return makeLayer('circle',recorded.circles,options);},
    circleMarker(_points,options){return makeLayer('circleMarker',recorded.circles,options);},
    marker(_latlng,options){return makeLayer('marker',recorded.markers,options);},
    divIcon(options){return options;},
    latLngBounds(points){return {points};},
  };
  const {PlannerMap}=await import('../../assets/js/ui/map.js');
  const project=createDemo();
  const map=new PlannerMap({innerHTML:''}, {onPick(){},onLegLength(){},onFinalLength(){},onHandle(){}});
  map.redraw(project,buildPlan(project));
  assert.equal(recorded.mapOptions.preferCanvas,false);
  assert.equal(recorded.renderers.length,1);
  assert.ok(recorded.paths.length>=3,`paths: ${recorded.paths.length}`);
  assert.ok(recorded.circles.length>=4,`circles: ${recorded.circles.length}`);
  for(const path of [...recorded.paths,...recorded.circles]){
    assert.equal(path.options.renderer,recorded.renderers[0]);
  }
  assert.ok(recorded.markers.length>=2);
  delete globalThis.L;
});

test('v0.1.3 selected leg is highlighted and a path click selects its leg',async()=>{
 const paths=[],circles=[],markers=[];let selected=null;
 const layer=(type,opts)=>({type,options:opts,handlers:{},addTo(){(type==='path'?paths:type==='circle'?circles:markers).push(this);return this;},on(k,fn){this.handlers[k]=fn;return this;},bindTooltip(content,options){this.tooltipContent=content;this.tooltipOptions=options;return this;}});
 globalThis.L={
  map(){return {setView(){return this;},on(){return this;},invalidateSize(){},getContainer(){return {style:{}};}};},
  svg(){return {svgRenderer:true};},tileLayer(){return {addTo(){}};},layerGroup(){return {addTo(){return this;},clearLayers(){}};},
  polyline(_pts,opts){return layer('path',opts);},circle(_pt,opts){return layer('circle',opts);},circleMarker(_pt,opts){return layer('circle',opts);},marker(_pt,opts){return layer('marker',opts);},divIcon(x){return x;}
 };
 const {PlannerMap}=await import('../../assets/js/ui/map.js');const project=previouslyPlannedProject();project.arrivalLegs=trafficPreset('left');
 project.display.referenceLineOpacity=.2;
 project.display.showCheckpointAltitudes=true;project.display.showDownwindMargin=true;
 const plan=buildPlan(project);
 const map=new PlannerMap({innerHTML:''},{onPick(){},onLegLength(){},onFinalLength(){},onHandle(){},onSelectLeg(id){selected=id;}});
 map.redraw(project,plan,project.departureLegs[0].id);
 assert.ok(paths.some(p=>p.options.color==='#f59e0b'&&p.options.weight>=5),'orange selected leg');
 assert.ok(paths.some(p=>p.options.color==='#fff'&&p.options.weight>=10),'selected leg halo');
 const base=paths.find(p=>p.options.color==='#2563eb');assert.ok(base?.handlers.click);
 base.handlers.click();assert.equal(selected,project.departureLegs[0].id);
 assert.equal(base.options.smoothFactor,0,'Leaflet simplification disabled for sampled arcs');
 assert.equal(circles.filter(c=>c.options.radius&&c.options.opacity===.2).length,8);
 assert.equal(paths.filter(p=>p.options.color==='#94a3b8'&&p.options.opacity===.2).length,8);
 assert.ok(circles.some(c=>String(c.tooltipContent).includes('高度余裕')));
 assert.ok(circles.some(c=>String(c.tooltipContent).includes('glp-cp-height')));
 assert.ok(circles.some(c=>/^<div class=\"glp-map-cp-content\">/.test(String(c.tooltipContent))));
 assert.ok(circles.every(c=>!String(c.tooltipContent).includes(' MSL</span>')),'permanent labels contain only values/units');
 delete globalThis.L;
});
test('v0.1.3 checkpoint altitudes and downwind margin labels can be disabled independently',async()=>{
 const circles=[];const marker=(opts)=>({options:opts,addTo(){circles.push(this);return this;},on(){return this;},bindTooltip(text,options){this.text=text;this.tooltipOptions=options;return this;}});
 globalThis.L={map(){return {setView(){return this;},on(){return this;},invalidateSize(){}};},svg(){return {};},tileLayer(){return {addTo(){}};},layerGroup(){return {addTo(){return this;},clearLayers(){}};},polyline(){return {addTo(){return this;},on(){return this;},bindTooltip(){return this;}};},circle(){return {addTo(){return this;}};},circleMarker(_pt,opts){return marker(opts);},marker(_pt,opts){return marker(opts);},divIcon(v){return v;}};
 const {PlannerMap}=await import('../../assets/js/ui/map.js');const project=createDemo();project.arrivalLegs=trafficPreset('right');project.display.showCheckpointAltitudes=false;project.display.showDownwindMargin=false;
 const map=new PlannerMap({innerHTML:''},{onPick(){},onLegLength(){},onFinalLength(){},onHandle(){}});
 map.redraw(project,buildPlan(project));
 assert.equal(circles.filter(c=>c.tooltipOptions?.permanent===true).length,0,'no permanent checkpoint labels when disabled');
 project.display.showDownwindMargin=true;circles.length=0;map.redraw(project,buildPlan(project));
 assert.equal(circles.filter(c=>c.tooltipOptions?.permanent===true).length,1,'only downwind label permanent');
 delete globalThis.L;
});
