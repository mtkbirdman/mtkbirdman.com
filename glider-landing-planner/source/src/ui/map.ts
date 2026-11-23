import type {ProjectV1,PlanResult,Point2,Side,Checkpoint,Sample} from '../model/types.js';
import {LocalFrame,headingVector,dot,sub,add,mul} from '../geometry/geo.js';
import {windVector} from '../flight/core.js';
declare const L:any;
type PickMode='departure'|'landing'|null;
interface Callbacks {
 onPick:(mode:'departure'|'landing',lat:number,lon:number)=>void;
 onLegLength:(side:Side,id:string,lengthM:number)=>void;
 onFinalLength:(lengthM:number)=>void;
 onHandle:(end:'start'|'end',lengthM:number)=>void;
 onSelectLeg?:(id:string)=>void;
}
const html=(value:unknown)=>String(value??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]??c));
const heightText=(m:number,unit:'m'|'ft')=>`${Math.round(unit==='ft'?m/0.3048:m).toLocaleString('ja-JP')} ${unit}`;
const marginText=(m:number,unit:'m'|'ft')=>`${m>0?'+':''}${heightText(m,unit)}`;

export class PlannerMap {
 private map:any;private graphics:any;private rings:any;private controls:any;private vectorRenderer:any;private picking:PickMode=null;private ready=false;
 constructor(private element:HTMLElement,private cb:Callbacks){
   if(typeof L==='undefined'){element.innerHTML='<div class="glp-map-unavailable">地図ライブラリの読込に失敗しました。インターネット接続・CDN制限をご確認ください。</div>';return;}
   this.map=L.map(element,{zoomControl:true,preferCanvas:false}).setView([35.4,136.7],14);
   // SVG remains necessary in Cocoon, and SVG simplification is disabled for smooth sampled arcs.
   this.vectorRenderer=L.svg({padding:0.5});
   L.tileLayer('https://cyberjapandata.gsi.go.jp/xyz/pale/{z}/{x}/{y}.png',{attribution:'&copy; <a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>',maxZoom:18,minZoom:3}).addTo(this.map);
   this.graphics=L.layerGroup().addTo(this.map);this.rings=L.layerGroup().addTo(this.map);this.controls=L.layerGroup().addTo(this.map);
   this.map.on('click',(ev:any)=>{if(this.picking){this.cb.onPick(this.picking,ev.latlng.lat,ev.latlng.lng);this.setPick(null);}});
   this.ready=true;setTimeout(()=>this.map.invalidateSize(),250);
 }
 setPick(mode:PickMode){this.picking=mode;if(this.ready)this.map.getContainer().style.cursor=mode?'crosshair':'';}
 fit(project:ProjectV1,plan:PlanResult){if(!this.ready)return;const frame=new LocalFrame(project.landing.point);const pts=[project.departure,project.landing.point,...[...plan.departure.samples,...(plan.connector?.samples??[]),...plan.arrival.samples,...plan.finalTarget].filter((_,i)=>i%Math.max(1,Math.ceil((plan.forwardAll.length+plan.requiredAll.length)/120))===0).map(p=>frame.toGeo(p))];this.map.fitBounds(L.latLngBounds(pts.map(p=>[p.latDeg,p.lonDeg])),{padding:[45,45],maxZoom:16});}
 redraw(project:ProjectV1,plan:PlanResult,selectedLegId:string|null=null){if(!this.ready)return;this.graphics.clearLayers();this.rings.clearLayers();this.controls.clearLayers();
 const frame=new LocalFrame(project.landing.point),loc=(p:Point2)=>{const g=frame.toGeo(p);return [g.latDeg,g.lonDeg]};
 const draw=(points:Point2[],color:string,weight=3,dashArray?:string,legId?:string)=>{
   if(points.length<2)return;
   const line=L.polyline(points.map(loc),{color,weight,opacity:.95,dashArray,lineCap:'round',lineJoin:'round',smoothFactor:0,renderer:this.vectorRenderer,bubblingMouseEvents:false});
   if(legId&&this.cb.onSelectLeg){line.on('click',()=>this.cb.onSelectLeg?.(legId));}
   line.addTo(this.graphics);
 };
 for(const seg of plan.departure.segments)draw(seg.samples,'#2563eb',3.5,undefined,seg.leg.id);
 if(plan.connector)draw(plan.connector.points,'#8b5cf6',3.5,plan.connector.feasible?undefined:'8,5');
 for(const seg of plan.arrival.segments)draw(seg.samples,'#059669',3.5,undefined,seg.leg.id);
 const end=plan.finalTarget[plan.finalTarget.length-1];
 if(end)draw([plan.finalStart,end],'#059669',4,'9,4');
 const selected=[...plan.departure.segments,...plan.arrival.segments].find(seg=>seg.leg.id===selectedLegId);
 if(selected&&selected.samples.length>1){
   // White halo + orange overlay remains visible on every basemap without changing the underlying leg colors.
   const pts=selected.samples.map(loc);
   L.polyline(pts,{color:'#fff',weight:11,opacity:.95,interactive:false,smoothFactor:0,lineJoin:'round',lineCap:'round',renderer:this.vectorRenderer}).addTo(this.graphics);
   L.polyline(pts,{color:'#f59e0b',weight:5.5,opacity:1,interactive:false,smoothFactor:0,lineJoin:'round',lineCap:'round',renderer:this.vectorRenderer}).addTo(this.graphics);
 }
 const depIcon=L.divIcon({className:'glp-map-point',html:'<span class="glp-point glp-point-blue">離</span>',iconSize:[25,25],iconAnchor:[12,12]});
 const landIcon=L.divIcon({className:'glp-map-point',html:'<span class="glp-point glp-point-green">着</span>',iconSize:[25,25],iconAnchor:[12,12]});
 L.marker([project.departure.latDeg,project.departure.lonDeg],{draggable:true,icon:depIcon,title:'離脱地点：ドラッグで変更'}).on('dragend',(e:any)=>{const x=e.target.getLatLng();this.cb.onPick('departure',x.lat,x.lng);}).addTo(this.controls);
 L.marker([project.landing.point.latDeg,project.landing.point.lonDeg],{draggable:true,icon:landIcon,title:'接地点：ドラッグで変更'}).on('dragend',(e:any)=>{const x=e.target.getLatLng();this.cb.onPick('landing',x.lat,x.lng);}).addTo(this.controls);
 const smallIcon=(name:string)=>L.divIcon({className:'glp-map-point',html:'<span class="glp-handle">'+html(name)+'</span>',iconSize:[19,19],iconAnchor:[9,9]});
 for(const side of ['departure','arrival'] as Side[]){const path=side==='departure'?plan.departure:plan.arrival;
 for(const segment of path.segments){if(segment.leg.type!=='straight')continue;const pt=side==='departure'?segment.end:segment.start;
 L.marker(loc(pt),{draggable:true,icon:smallIcon('↔'),title:'ドラッグで直線レグ長を変更（進行方向へ投影）'}).on('dragend',(ev:any)=>{
   const ll=ev.target.getLatLng(),grab=frame.toLocal({latDeg:ll.lat,lonDeg:ll.lng});
   const origin=side==='departure'?segment.start:segment.end;
   const unit=headingVector(segment.start.trackDeg);
   const vector=side==='departure'?sub(grab,origin):sub(origin,grab);
   const newLength=Math.max(1,Math.abs(dot(vector,unit)));
   this.cb.onLegLength(side,segment.leg.id,newLength);
 }).addTo(this.controls);
 }}
 L.marker(loc(plan.finalStart),{draggable:true,icon:smallIcon('F'),title:'ドラッグでファイナル長を変更'}).on('dragend',(ev:any)=>{const ll=ev.target.getLatLng(),p=frame.toLocal({latDeg:ll.lat,lonDeg:ll.lng});const direct=headingVector(project.landing.runwayTrackDeg);const v=sub(plan.landing,p);this.cb.onFinalLength(Math.max(20,dot(v,direct)));}).addTo(this.controls);
 if(plan.connector&&project.display.showHandles){
   const A=plan.departure.end,B=plan.arrival.start,da=headingVector(A.trackDeg),db=headingVector(B.trackDeg);
   const p1=add(A,mul(da,project.connector!.handleStartM)),p4=sub(B,mul(db,project.connector!.handleEndM));
   for(const [name,p,base,unit,end] of [['A',p1,A,da,'start'],['B',p4,B,mul(db,-1),'end']] as const){
      L.polyline([loc(base),loc(p)],{color:'#9333ea',dashArray:'3,5',weight:1.5,opacity:.8,smoothFactor:0,renderer:this.vectorRenderer}).addTo(this.controls);
      L.marker(loc(p),{draggable:true,icon:smallIcon(name),title:'Bezier制御ハンドル'}).on('dragend',(ev:any)=>{const ll=ev.target.getLatLng(),q=frame.toLocal({latDeg:ll.lat,lonDeg:ll.lng});this.cb.onHandle(end,Math.max(1,dot(sub(q,base),unit)));}).addTo(this.controls);
   }
 }
 // Merge checkpoints sharing one map coordinate, so downwind start is not covered by a redundant connector marker.
 const buckets=new Map<string,Checkpoint[]>();
 for(const cp of plan.checkpoints){const key=`${Math.round(cp.e/2)}:${Math.round(cp.n/2)}`;const existing=buckets.get(key)??[];existing.push(cp);buckets.set(key,existing);}
 let labelIndex=0;
 for(const cluster of buckets.values()){
   const cp=cluster.find(c=>c.checkpointKind==='downwind-start')??cluster[0];
   const downwind=cluster.find(c=>c.checkpointKind==='downwind-start');
   const showHeights=project.display.showCheckpointAltitudes;
   const showMargin=Boolean(downwind&&project.display.showDownwindMargin);
   const lines:string[]=[];
   if(showHeights){for(const c of cluster){const kind=c.altitudeType==='predicted'?'予想高度':c.altitudeType==='required'?'必要高度':'目標高度';lines.push(`<span>${html(c.label)}：${kind} ${html(heightText(c.altitudeM,project.display.altitudeUnit))} MSL</span>`);}}
   if(showMargin&&downwind){
     const value=downwind.marginM===undefined?'評価不可':html(marginText(downwind.marginM,project.display.altitudeUnit));
     lines.push(`<strong>ダウンウィンド開始｜高度余裕 ${value}</strong>`);
   }
   const hover=cluster.map(c=>`${c.label}: ${heightText(c.altitudeM,project.display.altitudeUnit)} MSL`).join(' / ');
   const marker=L.circleMarker(loc(cp),{radius:4,fillColor:downwind?'#b45309':'#d97706',fillOpacity:1,color:'#fff',weight:1,interactive:true,bubblingMouseEvents:false,renderer:this.vectorRenderer});
   marker.bindTooltip(lines.length?`<div class="glp-map-cp-content">${lines.join('')}</div>`:html(hover),
     lines.length?{permanent:true,direction:labelIndex%2?'bottom':'top',offset:[0,labelIndex%2?8:-5],opacity:0.96,className:'glp-checkpoint-label',interactive:false}:{permanent:false});
   marker.addTo(this.graphics);labelIndex++;
 }
 if(project.display.showDistanceRings){for(const radius of [250,500,1000,2000])L.circle([project.landing.point.latDeg,project.landing.point.lonDeg],{radius,color:'#7a8797',weight:1,opacity:project.display.distanceRingOpacity,fill:false,dashArray:'4,5',interactive:false,renderer:this.vectorRenderer}).addTo(this.rings);}
 if(project.display.showBearingLines){for(let i=0;i<360;i+=45){const vec=headingVector(i),a=loc(plan.landing),b=loc(add(plan.landing,mul(vec,2000)));L.polyline([a,b],{color:'#94a3b8',weight:.8,opacity:project.display.bearingLineOpacity,dashArray:'3,7',smoothFactor:0,interactive:false,renderer:this.vectorRenderer}).addTo(this.rings);}}
 if(project.display.showWindArrow&&project.wind.speedMps>0){const v=windVector(project.wind),unit=mul(v,1/Math.hypot(v.e,v.n)),center=add(plan.landing,{e:150,n:150}),tail=add(center,mul(unit,-80)),head=add(center,mul(unit,80));L.polyline([loc(tail),loc(head)],{color:'#c2410c',weight:3,renderer:this.vectorRenderer}).bindTooltip('風: '+project.wind.fromDeg+'°から '+project.wind.speedMps+' m/s').addTo(this.rings);}
 }
}
