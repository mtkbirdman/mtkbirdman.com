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

test('wind arrow points downwind in all eight directions and obeys display controls',async()=>{
 const {LocalFrame}=await import('../../assets/js/geometry/geo.js');
 const {PlannerMap}=await import('../../assets/js/ui/map.js');
 const paths=[];const prev=globalThis.L;
 const layer=(points,options)=>({points,options,addTo(){if(options.color==='#c2410c')paths.push(this);return this;},on(){return this;},bindTooltip(text){this.tooltip=text;return this;}});
 globalThis.L={map(){return {setView(){return this;},on(){return this;},invalidateSize(){}};},svg(){return {};},tileLayer(){return {addTo(){}};},layerGroup(){return {addTo(){return this;},clearLayers(){}};},polyline:layer,circle:layer,circleMarker:layer,marker:layer,divIcon(x){return x;}};
 try{
  const p=createDemo(),frame=new LocalFrame(p.landing.point);
  const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(){}});
  const local=([latDeg,lonDeg])=>frame.toLocal({latDeg,lonDeg});
  const near=(a,b)=>assert.ok(Math.abs(a-b)<.002,`${a} != ${b}`);
  for(let fromDeg=0;fromDeg<360;fromDeg+=45){
   p.wind={fromDeg,speedMps:10};p.display.magneticDeclinationDeg=15;paths.length=0;
   map.redraw(p,buildPlan(p));assert.equal(paths.length,2);
   const [shaft,head]=paths;assert.equal(shaft.points.length,2);assert.equal(head.points.length,3);
   assert.deepEqual(head.points[1],shaft.points[1]);
   assert.match(shaft.tooltip,new RegExp(`${fromDeg}°から 10 m/s`));
   const [a,b]=shaft.points.map(local),[l,t,r]=head.points.map(local);
   const e=-Math.sin(fromDeg*Math.PI/180),n=-Math.cos(fromDeg*Math.PI/180);
   near(b.e-a.e,160*e);near(b.n-a.n,160*n);
   for(const wing of [l,r])near((wing.e-t.e)*e+(wing.n-t.n)*n,-45);
   near((l.e-t.e)*(-n)+(l.n-t.n)*e,32);
   near((r.e-t.e)*(-n)+(r.n-t.n)*e,-32);
   assert.equal(head.options.renderer,shaft.options.renderer);
  }
  paths.length=0;p.wind.speedMps=0;map.redraw(p,buildPlan(p));assert.equal(paths.length,0);
  p.wind.speedMps=10;p.display.showWindArrow=false;map.redraw(p,buildPlan(p));assert.equal(paths.length,0);
 }finally{globalThis.L=prev;}
});

test('GLP-015/016/017: invalid preview, partial selected routes, fit and bounded drag',async()=>{
 const {makeStraight}=await import('../../assets/js/model/defaults.js');
 const {LocalFrame}=await import('../../assets/js/geometry/geo.js');
 const {PlannerMap}=await import('../../assets/js/ui/map.js');
 const prevL=globalThis.L,prevW=globalThis.window;
 const layers=[],fits=[],calls=[];
 const make=(kind)=>(points,options)=>({kind,points,options,handlers:{},addTo(){layers.push(this);return this},on(k,f){this.handlers[k]=f;return this},bindTooltip(){return this}});
 globalThis.window={matchMedia:()=>({matches:true})};
 globalThis.L={map(){return {setView(){return this},on(){return this},invalidateSize(){},fitBounds(v){fits.push(v)}}},svg(){return {}},tileLayer(){return {addTo(){}}},layerGroup(){return {addTo(){return this},clearLayers(){}}},polyline:make('line'),marker:make('marker'),circle:make('circle'),circleMarker:make('checkpoint'),divIcon(v){return v},latLngBounds(v){return v}};
 try{
   const p=createDemo();p.landing.point={latDeg:0,lonDeg:0};p.landing.runwayTrackDeg=90;p.arrivalLegs=[makeStraight(300)];p.departure={...p.departure,latDeg:0,lonDeg:Math.asin(-1900/6378137)*180/Math.PI,headingDeg:90};p.connector={type:'bezier_connector',handleStartM:200,handleEndM:590.6251,curvaturePolicy:'match-neighbors',adjustmentMode:'manual'};p.display.showHandles=true;
   const plan=buildPlan(p);assert.equal(plan.connector.samples.length,0);const map=new PlannerMap({innerHTML:''},{onPick(){},onHandle(...a){calls.push(a)},onSelectLeg(id){calls.push(id)}});
   map.redraw(p,plan);map.fit(p,plan);assert.ok(layers.some(l=>l.kind==='line'&&l.options.color==='#8b5cf6'&&l.options.dashArray==='8,5'));assert.ok(fits.at(-1).length>4);
   const handles=layers.filter(l=>l.kind==='marker'&&l.options.title==='Bezier制御ハンドル');assert.equal(handles.length,2);
   const frame=new LocalFrame(p.landing.point),A=plan.departure.end;
   for(const e of [A.e-100,A.e+55000,NaN]){const q=Number.isFinite(e)?frame.toGeo({e,n:A.n}):{latDeg:NaN,lonDeg:NaN};handles[0].handlers.dragend({target:{getLatLng(){return {lat:q.latDeg,lng:q.lonDeg}}}});}
   assert.deepEqual(calls.slice(0,2),[['start',1],['start',50000]]);assert.equal(calls[2][0],'start');assert.ok(!Number.isFinite(calls[2][1]));
   p.connector=null;p.arrivalLegs=[];p.departure={...p.departure,latDeg:.442,lonDeg:0,headingDeg:0};p.departureLegs=[makeStraight(100),makeStraight(2000)];const partial=buildPlan(p);const id=p.departureLegs[1].id;layers.length=0;map.redraw(p,partial,id);
   for(const color of ['#2563eb','#fff','#f59e0b'])assert.ok(layers.some(l=>l.kind==='line'&&l.options.color===color&&l.options.dashArray==='8,5'));
   assert.ok(layers.some(l=>l.kind==='line'&&l.options.color==='#2563eb'&&!l.options.dashArray&&l.options.weight===3.5));
   const hit=layers.filter(l=>l.options.weight===22).at(-1);assert.ok(hit);hit.handlers.click({});assert.equal(calls.at(-1),id);
 }finally{globalThis.L=prevL;globalThis.window=prevW;}
});
