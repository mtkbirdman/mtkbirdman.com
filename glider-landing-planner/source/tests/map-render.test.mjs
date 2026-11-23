import test from 'node:test';
import assert from 'node:assert/strict';
import {createDemo} from '../../assets/js/model/defaults.js';
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
