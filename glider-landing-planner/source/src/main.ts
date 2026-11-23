import type {ProjectV1,FlightLeg,Side,PlanResult} from './model/types.js';
import {createDemo,makeStraight,makeTurn,trafficPreset,ENGINE_VERSION} from './model/defaults.js';
import {buildPlan} from './planning/plan.js';
import {chooseConnector} from './planning/connector.js';
import {validateProject,parseProject,downloadProject} from './io/project.js';
import {saveSnapshot,latestSnapshot} from './io/autosave.js';
import {PlannerMap} from './ui/map.js';
import {chartSvg,fmtAlt,altitudeTimeSeries,fmtElapsed} from './ui/chart.js';
import {openPrintReport} from './ui/print.js';
import {checkpointTable} from './ui/checkpoints.js';
import {heightInputText,FEET_PER_METER} from './ui/height-units.js';
const escapeHTML=(val:unknown)=>String(val??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]??c));
const n=(value:number,d=2)=>Number.isFinite(value)?Number(value.toFixed(d)):0;
const label=(s:string)=>escapeHTML(s);
const options=(items:[string,string][],selected:string)=>items.map(([v,t])=>`<option value="${v}" ${v===selected?'selected':''}>${t}</option>`).join('');
const engineNote='教育・机上検討用。実際の飛行可否、地形・障害物、フレア・スポイラー性能は評価しません。';
class App {
 private project:ProjectV1=createDemo();private plan:PlanResult=buildPlan(this.project);
 private history:ProjectV1[]=[];private future:ProjectV1[]=[];private map:PlannerMap;
 private selected:string|null=null;private pick:'departure'|'landing'|null=null;private saveTimer:number|undefined;
 private status='未保存';private notice='';private editRevision=0;
 constructor(private root:HTMLElement){
   root.innerHTML=`<div class="glp-shell"><div class="glp-toolbar"></div><div class="glp-app-body"><aside class="glp-sidebar"><div class="glp-sidebar-scroll"></div></aside><section class="glp-main"><div class="glp-mapbar"><div class="glp-map-actions"><button data-action="pick-departure" class="glp-btn glp-btn-quiet">離脱地点を指定</button><button data-action="pick-landing" class="glp-btn glp-btn-quiet">着陸地点を指定</button><button data-action="fit" class="glp-btn glp-btn-quiet">全体表示</button></div><span class="glp-map-hint">地点・直線端点・紫のハンドルはドラッグできます</span></div><div class="glp-map"></div><div class="glp-overview"><div class="glp-overview-top"><strong>高度プロファイル</strong><span class="glp-chart-legend"><i class="glp-line-blue"></i>離脱側予想高度 <i class="glp-line-green"></i>着陸側逆算必要高度</span></div><div class="glp-chart"></div><div class="glp-metrics"></div><div class="glp-checkpoints"></div><div class="glp-diagnostics"></div></div></section></div><div class="glp-footnote">${engineNote}　｜　地図：<a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>（地理院タイル）</div><input class="glp-import-file" type="file" accept=".json,.glp.json,application/json" hidden></div>`;
   this.map=new PlannerMap(root.querySelector('.glp-map') as HTMLElement,{
     onPick:(type,lat,lon)=>{this.pick=null;this.map.setPick(null);this.edit(p=>{if(type==='departure'){p.departure.latDeg=lat;p.departure.lonDeg=lon;}else{p.landing.point.latDeg=lat;p.landing.point.lonDeg=lon;}});},
     onLegLength:(side,id,lengthM)=>this.edit(p=>{const l=(side==='departure'?p.departureLegs:p.arrivalLegs).find(v=>v.id===id);if(l?.type==='straight')l.lengthM=lengthM;}),
     onFinalLength:lengthM=>this.edit(p=>{p.landing.final.groundLengthM=lengthM;}),
     onHandle:(end,lengthM)=>this.edit(p=>{if(p.connector){p.connector.adjustmentMode='manual';if(end==='start')p.connector.handleStartM=lengthM;else p.connector.handleEndM=lengthM;}}),
     onSelectLeg:(id)=>{this.selected=this.selected===id?null:id;this.renderSidebar();this.map.redraw(this.project,this.plan,this.selected);this.updateCheckpointSelection();}
   });
   root.addEventListener('click',e=>{const btn=(e.target as HTMLElement).closest<HTMLButtonElement>('button[data-action]');if(btn){e.preventDefault();this.handleAction(btn.dataset.action??'',btn.dataset);}});
   root.addEventListener('change',e=>{const el=e.target as HTMLInputElement|HTMLSelectElement;
      if(el.matches('[data-path]'))this.handleChange(el);
      else if(el.classList.contains('glp-import-file'))void this.doImport(el as HTMLInputElement);
   });
   this.render();this.map.fit(this.project,this.plan);
   const revisionAtStartup=this.editRevision;
   void latestSnapshot().then(p=>{
     if(!p||this.editRevision!==revisionAtStartup)return; // Don't overwrite edits made while IndexedDB was loading.
     let valid:ProjectV1;
     try{valid=validateProject(p);}catch(e){
       this.notice='自動保存データが不正なため復元できません: '+(e instanceof Error?e.message:String(e));
       this.renderResults();return;
     }
     this.project=valid;this.plan=buildPlan(valid);this.history=[];this.future=[];this.selected=null;
     this.status='自動保存から復元しました';this.render();this.map.fit(this.project,this.plan);
     // Persist migrated old bank angles once; failure must not discard the restored plan.
     void saveSnapshot(this.project).catch(()=>{this.status='自動保存に失敗（JSON保存を推奨）';this.renderToolbar();});
   }).catch(()=>{ /* IndexedDB may be blocked or contain invalid legacy data: keep demo visible. */ });
 }
 private getField(path:string):unknown{let a:any=this.project;for(const k of path.split('.'))a=a[k];return a;}
 private control(path:string,title:string,unit='',step='any',scale=1):string {
   const raw=this.getField(path);
   const isHeight=path==='departure.altitudeMSLM'||path==='landing.elevationMSLM';
   const heightUnit=this.project.display.altitudeUnit;
   const actualScale=isHeight&&heightUnit==='ft'?FEET_PER_METER:scale;
   const fieldUnit=isHeight?`${heightUnit}, MSL`:unit;
   const value=isHeight&&typeof raw==='number'?heightInputText(raw,heightUnit):typeof raw==='number'?n(raw*scale,scale===3.6?2:5):raw??'';
   return `<label class="glp-field"><span>${title}${fieldUnit?` <em>${fieldUnit}</em>`:''}</span><input type="number" step="${step}" value="${escapeHTML(value)}" data-path="${path}" data-scale="${actualScale}"></label>`;
 }
 private checkbox(path:string,title:string):string{return `<label class="glp-check"><input type="checkbox" data-path="${path}" ${this.getField(path)?'checked':''}><span>${title}</span></label>`;}
 private opacityControl(path:string,title:string):string{const value=Math.round(Number(this.getField(path))*100);return `<label class="glp-opacity-field"><span>${title} <output>${value}%</output></span><input type="range" min="0" max="100" step="5" value="${value}" data-path="${path}" data-scale="100"></label>`;}
 private select(path:string,title:string,items:[string,string][]):string{return `<label class="glp-field"><span>${title}</span><select data-path="${path}">${options(items,String(this.getField(path)))}</select></label>`;}
 private legList(side:Side){const legs=side==='departure'?this.project.departureLegs:this.project.arrivalLegs;
   const modeOpts:[string,string][]=[['heading','機首方位指定（風で偏流）'],['ground-track','対地航跡指定（風を補正）']];
   const lines=legs.map((leg,i)=>{
     const active=this.selected===leg.id;const title=leg.type==='straight'?`直線 ${Math.round(leg.lengthM)} m`:`${leg.direction==='right'?'右':'左'}旋回 ${leg.angleDeg}° / バンク${leg.bankDeg}°`;
     const group=side==='departure'?'departureLegs':'arrivalLegs';
     const base=`${group}.${i}`;
     return `<div class="glp-leg ${active?'glp-leg-active':''}"><div class="glp-leg-line"><button class="glp-leg-title" data-action="select-leg" data-id="${leg.id}" title="編集"><span class="glp-leg-symbol">${leg.type==='straight'?'━':'⤷'}</span><span>${i+1}. ${label(title)}</span></button><span class="glp-leg-buttons"><button title="上に移動" data-action="move-up" data-side="${side}" data-id="${leg.id}">↑</button><button title="下に移動" data-action="move-down" data-side="${side}" data-id="${leg.id}">↓</button><button title="削除" data-action="delete-leg" data-side="${side}" data-id="${leg.id}">×</button></span></div>
     ${active?`<div class="glp-leg-edit"><div class="glp-two">${this.select(base+'.guidance','風補正方式',modeOpts)}${this.control(base+'.speedOverrideMps','個別速度（空欄＝共通）','km/h','0.1',3.6)}</div>
       ${leg.type==='straight'?`<div class="glp-two">${this.control(base+'.lengthM','直線距離','m')}${this.control(base+'.directionDeg','指定方位（空欄＝継承）','°')}</div>`:`<div class="glp-two">${this.select(base+'.direction','旋回方向',[['right','右'],['left','左']])}${this.control(base+'.angleDeg','旋回角（多周旋回可）','°')}</div><div class="glp-two">${this.control(base+'.bankDeg','バンク角','°')}</div>`}
       <label class="glp-field"><span>チェックポイント名</span><input type="text" maxlength="70" value="${escapeHTML(leg.checkpointLabel??'')}" data-path="${base}.checkpointLabel"></label></div>`:''}</div>`;
   }).join('');
   const adds=`<div class="glp-adds"><button class="glp-btn glp-btn-sub" data-action="add-straight" data-side="${side}">＋ 直線</button><button class="glp-btn glp-btn-sub" data-action="add-turn" data-side="${side}">＋ 旋回</button></div>`;
   const list=`<div class="glp-leg-list">${lines||'<p class="glp-muted">まだレグがありません。</p>'}</div>`;
   return side==='arrival'?adds+list:list+adds;
 }
 private renderSidebar(){
  const sidebarScroll=this.root.querySelector('.glp-sidebar-scroll') as HTMLElement;
  const preservedScroll=sidebarScroll.scrollTop;
  sidebarScroll.innerHTML=`
   <div class="glp-panel-heading"><div class="glp-eyebrow">FLIGHT PARAMETERS</div><h2>飛行条件・経路設計</h2></div>
   <div class="glp-panel"><h3>基本条件</h3><div class="glp-two">${this.control('defaultSpeedMps','共通対気速度','km/h','0.1',3.6)}${this.control('maxBankDeg','最大バンク','°')}</div>
   ${this.select('performance.model','沈下計算モデル',[['fixed-ld','A：固定L/D'],['quadratic-drag','B：簡易抗力モデル']])}
   <div class="glp-two">${this.control('performance.fixedLD','モデルA：L/D')}${this.control('performance.bestLD','モデルB：最大L/D')}</div>${this.control('performance.bestLDAtMps','モデルB：最大L/D速度','km/h','0.1',3.6)}
   <div class="glp-two">${this.control('wind.fromDeg','風向','°')}${this.control('wind.speedMps','風速','m/s')}</div>
   ${this.select('display.altitudeUnit','高度入力・表示単位（MSL）',[['ft','ft'],['m','m']])}</div>
   <div class="glp-panel"><h3><span class="glp-dot blue"></span>離脱地点・条件</h3><div class="glp-two">${this.control('departure.latDeg','緯度','°')}${this.control('departure.lonDeg','経度','°')}</div><div class="glp-two">${this.control('departure.altitudeMSLM','離脱高度')}${this.control('departure.headingDeg','開始機首方位','°')}</div><p class="glp-tip">地図の「離脱地点を指定」または青いマーカーのドラッグで変更できます。</p></div>
   <div class="glp-panel"><h3><span class="glp-dot green"></span>着陸地点・ファイナル</h3><div class="glp-two">${this.control('landing.point.latDeg','緯度','°')}${this.control('landing.point.lonDeg','経度','°')}</div><div class="glp-two">${this.control('landing.elevationMSLM','着陸標高')}${this.control('landing.runwayTrackDeg','真方位（進入）','°')}</div>${this.control('landing.final.pathAngleDeg','目標パス角','°')}<p class="glp-tip">最終パス角は目標値であり、スポイラー・フレアは未モデル化です。</p></div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot blue"></span>離脱側レグ（順方向）</h3><small>${this.project.departureLegs.length}件</small></div>${this.legList('departure')}</div>
   <div class="glp-panel"><div class="glp-flex-between"><h3><span class="glp-dot green"></span>着陸側レグ（逆方向編集）</h3><small>${this.project.arrivalLegs.length}件</small></div>
   <div class="glp-preset"><button class="glp-btn glp-btn-soft" data-action="preset-left">左場周プリセット</button><button class="glp-btn glp-btn-soft" data-action="preset-right">右場周プリセット</button></div>${this.control('landing.final.groundLengthM','直線ファイナル長','m')}${this.legList('arrival')}
   <p class="glp-tip">新しいレグは最上流へ追加され、一覧は実際の飛行順です。</p></div>
   <div class="glp-panel"><h3><span class="glp-dot purple"></span>滑らか接続（5次Bezier）</h3>
   ${this.checkbox('display.showHandles','地図に制御ハンドルを表示')}
   ${this.project.connector?this.select('connector.adjustmentMode','接続曲線の調整',[['auto','自動調整（標準）'],['manual','手動調整（ハンドル操作）']]):'<p class="glp-muted">接続がオフになっています。</p>'}
   <div class="glp-adds"><button class="glp-btn glp-btn-primary" data-action="optimize">曲線を自動調整</button><button class="glp-btn glp-btn-sub" data-action="toggle-connector">${this.project.connector?'接続を解除':'接続を作成'}</button></div><p class="glp-tip">接線と曲率を一致させます。大きいバンク角や風で不成立になる場合があります。</p></div>
   <div class="glp-panel"><h3>地図・表示</h3>${this.checkbox('display.showDistanceRings','着陸点から距離リング')}${this.checkbox('display.showBearingLines','方位線')}${this.checkbox('display.showWindArrow','風向矢印')}${this.checkbox('display.showDownwindMargin','ダウンウィンド開始の高度余裕を表示')}${this.checkbox('display.showCheckpointAltitudes','各チェックポイント高度を地図に表示')}${this.opacityControl('display.referenceLineOpacity','補助線の不透明度')}<div class="glp-two">${this.control('display.distanceRingIntervalM','距離リング間隔','m','25')}${this.select('display.bearingLineIntervalDeg','方位線間隔',[['15','15°'],['30','30°'],['45','45°'],['60','60°'],['90','90°'],['120','120°'],['180','180°']])}</div>${this.control('display.magneticDeclinationDeg','磁気偏角（東偏＋／西偏−）','°','0.1')}<p class="glp-tip">磁気偏角は補助方位線のみへ反映。飛行計算・風向・滑走路方位は真方位基準です。</p></div>`;
  sidebarScroll.scrollTop=preservedScroll;
 }
 private renderToolbar(){(this.root.querySelector('.glp-toolbar') as HTMLElement).innerHTML=`<div class="glp-brand"><div class="glp-logo">GLP</div><div><strong>Glider Landing Planner</strong><span>飛行経路シミュレーター / テスト版 v0.1.6</span></div></div><div class="glp-tool-group">
  <input type="text" class="glp-project-name" data-path="name" value="${escapeHTML(this.project.name)}" aria-label="プロジェクト名" maxlength="120">
  <button class="glp-btn glp-btn-quiet" data-action="new">新規</button><button class="glp-btn glp-btn-quiet" data-action="undo" ${!this.history.length?'disabled':''}>元に戻す</button><button class="glp-btn glp-btn-quiet" data-action="redo" ${!this.future.length?'disabled':''}>やり直す</button><button class="glp-btn glp-btn-quiet" data-action="export">JSON保存</button><button class="glp-btn glp-btn-quiet" data-action="import">JSON読込</button><button class="glp-btn glp-btn-primary" data-action="print">PDF印刷</button></div><div class="glp-save-state">● ${escapeHTML(this.status)}</div>`;}
 private renderResults(){const p=this.plan,unit=this.project.display.altitudeUnit;
 const downwind=p.checkpoints.find(c=>c.checkpointKind==='downwind-start');
 const connected=p.connector?.feasible&&downwind?.marginM!==undefined;
 const totalKm=(p.forwardAll.at(-1)?.sM??0)/1000;
 const gapM=p.connector?p.departure.end.h-p.connector.lossM-p.arrival.start.h:null;
 const flightTime=altitudeTimeSeries(p);
 const m=`<div class="glp-stat"><small>離脱高度 MSL</small><strong>${fmtAlt(this.project.departure.altitudeMSLM,unit)}</strong></div><div class="glp-stat"><small>推定総距離</small><strong>${Number.isFinite(totalKm)?totalKm.toFixed(2):'—'} km</strong></div><div class="glp-stat"><small>推定飛行時間${flightTime.timingProvisional?'（暫定）':''}</small><strong>${fmtElapsed(flightTime.durationS)}</strong></div><div class="glp-stat"><small>高度余裕</small><strong>${p.connector?.feasible&&gapM!==null?fmtAlt(gapM,unit):'評価不可'}</strong></div>`;
 (this.root.querySelector('.glp-chart') as HTMLElement).innerHTML=chartSvg(p,this.project);
 (this.root.querySelector('.glp-metrics') as HTMLElement).innerHTML=m;
 (this.root.querySelector('.glp-checkpoints') as HTMLElement).innerHTML=`<div class="glp-overview-top"><strong>チェックポイント高度</strong><small>滑空場基準高度（着陸標高を0とする相対高度）／MSL高度</small></div><div class="glp-tablewrap">${checkpointTable(p,this.project)}</div><p class="glp-tip">高度余裕は接続点での予想・必要高度の差です。着陸地点基準の高さではありません。</p>${downwind&&!connected?'<p class="glp-tip">ダウンウィンド開始の高度余裕は、接続条件が成立していないため評価不可です。</p>':''}`;
 const distinct=new Map<string,string>();for(const notice of p.notices)distinct.set(notice.code+notice.message,notice.message);
 const top=p.notices.filter(i=>i.level==='error');const warnings=p.notices.filter(i=>i.level==='warning');
 (this.root.querySelector('.glp-diagnostics') as HTMLElement).innerHTML=`<div class="glp-overview-top"><strong>計算上の診断</strong><span class="glp-status-chip ${top.length?'glp-status-bad':warnings.length?'glp-status-warn':'glp-status-note'}">${top.length?'不成立の条件あり':warnings.length?'警告あり':'近似モデル'}</span></div>
 <div class="glp-diagnostic-bank">Bezier接続最大必要バンク：<strong>${p.connector?p.connector.maxBankDeg.toFixed(1)+'°':'—'}</strong>（設定上限 ${this.project.maxBankDeg.toFixed(1)}°）</div>${this.notice?`<div class="glp-inline-message">${label(this.notice)}</div>`:''}<div class="glp-alerts">${[...distinct.values()].slice(0,12).map(t=>`<div class="glp-alert">${label(t)}</div>`).join('')||'<p>計算エラーは報告されていません。</p>'}</div>`;
 this.updateCheckpointSelection();
 }
 private updateCheckpointSelection(){
   for(const row of this.root.querySelectorAll<HTMLElement>('.glp-checkpoints tr[data-leg-id]')){
     row.classList.toggle('glp-checkpoint-active',!!this.selected&&row.dataset.legId===this.selected);
   }
 }
 private printReport(){
   const mapEl=this.root.querySelector<HTMLElement>('.glp-map');
   if(!mapEl)throw new Error('地図が見つかりません');
   openPrintReport({project:this.project,plan:this.plan,mapEl,
     metricsHtml:(this.root.querySelector('.glp-metrics') as HTMLElement).innerHTML,
     diagnosticsHtml:(()=>{const diag=this.root.querySelector('.glp-diagnostics') as HTMLElement;return ['.glp-diagnostic-bank','.glp-inline-message','.glp-alerts'].map(sel=>diag.querySelector(sel)?.outerHTML??'').join('');})(),
     title:'Glider Landing Planner',engineVersion:ENGINE_VERSION});
 }
 private render(){this.renderToolbar();this.renderSidebar();this.renderResults();this.map.redraw(this.project,this.plan,this.selected);
 const bb=this.root.querySelectorAll<HTMLButtonElement>('.glp-map-actions button');bb.forEach(b=>b.classList.toggle('glp-active',b.dataset.action==='pick-'+this.pick));}
 private optimizationSignature(p:ProjectV1):string{
   return JSON.stringify({speed:p.defaultSpeedMps,bank:p.maxBankDeg,performance:p.performance,wind:p.wind,
     departure:p.departure,landing:p.landing,departureLegs:p.departureLegs,arrivalLegs:p.arrivalLegs});
 }
 private recalculate(p:ProjectV1,previous:ProjectV1):PlanResult{
   const preliminary=buildPlan(p);
   if(!p.connector||p.connector.adjustmentMode!=='auto')return preliminary;
   const modeChanged=previous.connector?.adjustmentMode!=='auto';
   const shapeChanged=this.optimizationSignature(previous)!==this.optimizationSignature(p);
   if(!modeChanged&&!shapeChanged)return preliminary;
   const optimized=chooseConnector(preliminary.departure.end,preliminary.arrivalStart,p);
   p.connector.handleStartM=optimized.handleStartM;
   p.connector.handleEndM=optimized.handleEndM;
   return buildPlan(p);
 }
 private edit(mutate:(p:ProjectV1)=>void){this.editRevision++;const p=structuredClone(this.project),original=structuredClone(this.project);
 try{mutate(p);validateProject(p);const next=this.recalculate(p,original);this.history.push(original);if(this.history.length>50)this.history.shift();this.future=[];this.project=p;this.plan=next;this.notice='';this.status='未保存';this.render();this.scheduleSave();}
 catch(e){this.notice='変更を適用できません: '+(e instanceof Error?e.message:String(e));this.renderResults();}}
 private scheduleSave(){if(this.saveTimer)window.clearTimeout(this.saveTimer);this.saveTimer=window.setTimeout(()=>{this.status='保存中…';this.renderToolbar();void saveSnapshot(this.project).then(()=>{this.status='自動保存済';this.renderToolbar();}).catch(()=>{this.status='自動保存失敗（JSONを保存してください）';this.renderToolbar();});},1000);}
 private handleChange(input:HTMLInputElement|HTMLSelectElement){const path=input.dataset.path;if(!path)return;
 let value:any;
 if(input instanceof HTMLInputElement&&input.type==='checkbox')value=input.checked;
 else if(input instanceof HTMLInputElement&&(input.type==='number'||input.type==='range')){
   if(input.value.trim()===''){if(path.endsWith('speedOverrideMps')||path.endsWith('directionDeg'))value=null;else {this.notice='数値項目を空欄にはできません。';this.render();return;}}
   else{value=Number(input.value)/Number(input.dataset.scale??1);if(!Number.isFinite(value)){this.notice='数値が不正です。';this.render();return;}}
 }else if(input instanceof HTMLSelectElement&&path==='display.bearingLineIntervalDeg')value=Number(input.value);
 else value=input.value;
 this.edit(p=>{const keys=path.split('.');let target:any=p;for(let i=0;i<keys.length-1;i++)target=target[keys[i]];if(path.endsWith('checkpointLabel')&&!value)delete target[keys.at(-1)!];else target[keys.at(-1)!]=value;});
 }
 private handleAction(action:string,data:DOMStringMap){if(action==='fit'){this.map.fit(this.project,this.plan);return;}
 if(action==='pick-departure'||action==='pick-landing'){this.pick=action==='pick-departure'?'departure':'landing';this.map.setPick(this.pick);this.render();return;}
 if(action==='import'){(this.root.querySelector('.glp-import-file') as HTMLInputElement).click();return;}
 if(action==='export'){try{downloadProject(this.project);this.status='JSONを書き出しました';this.renderToolbar();}catch(e){this.notice=String(e);this.renderResults();}return;}
 if(action==='print'){try{this.printReport();}catch(e){this.notice='PDF印刷を開始できません: '+(e instanceof Error?e.message:String(e));this.renderResults();}return;}
 if(action==='undo'){if(this.history.length){this.editRevision++;this.future.push(structuredClone(this.project));this.project=this.history.pop()!;this.plan=buildPlan(this.project);this.render();this.scheduleSave();}return;}
 if(action==='redo'){if(this.future.length){this.editRevision++;this.history.push(structuredClone(this.project));this.project=this.future.pop()!;this.plan=buildPlan(this.project);this.render();this.scheduleSave();}return;}
 if(action==='new'){if(window.confirm('現在の計画を破棄して新規計画を作成しますか？（先にJSON保存を推奨）')){this.edit(p=>Object.assign(p,createDemo()));this.map.fit(this.project,this.plan);}return;}
 if(action==='select-leg'){this.selected=this.selected===data.id?null:data.id??null;this.renderSidebar();this.map.redraw(this.project,this.plan,this.selected);this.updateCheckpointSelection();return;}
 if(action==='add-straight'||action==='add-turn'){
 const leg=action==='add-straight'?makeStraight():makeTurn();const side=data.side as Side;
 this.edit(p=>{if(side==='departure')p.departureLegs.push(leg);else p.arrivalLegs.unshift(leg);});this.selected=leg.id;this.renderSidebar();this.map.redraw(this.project,this.plan,this.selected);this.updateCheckpointSelection();return;}
 if(action==='delete-leg'||action==='move-up'||action==='move-down'){const side=data.side as Side,id=data.id;
 this.edit(p=>{const list=side==='departure'?p.departureLegs:p.arrivalLegs;const i=list.findIndex(l=>l.id===id);if(i<0)return;
 if(action==='delete-leg')list.splice(i,1);
 else {const dest=i+(action==='move-up'?-1:1);if(dest<0||dest>=list.length)return;[list[i],list[dest]]=[list[dest],list[i]];}});
 if(action==='delete-leg'&&this.selected===id)this.selected=null;this.renderSidebar();this.map.redraw(this.project,this.plan,this.selected);this.updateCheckpointSelection();return;}
 if(action==='preset-left'||action==='preset-right'){
 if(this.project.arrivalLegs.length&&!window.confirm('着陸側レグをプリセットで置き換えますか？'))return;
 this.edit(p=>{p.arrivalLegs=trafficPreset(action==='preset-left'?'left':'right');});return;}
 if(action==='toggle-connector'){this.edit(p=>{p.connector=p.connector?null:{type:'bezier_connector',handleStartM:Math.max(20,Math.hypot(this.plan.arrivalStart.e-this.plan.departure.end.e,this.plan.arrivalStart.n-this.plan.departure.end.n)*.35),handleEndM:Math.max(20,Math.hypot(this.plan.arrivalStart.e-this.plan.departure.end.e,this.plan.arrivalStart.n-this.plan.departure.end.n)*.35),curvaturePolicy:'match-neighbors',adjustmentMode:'auto'};});return;}
 if(action==='optimize'){const result=chooseConnector(this.plan.departure.end,this.plan.arrivalStart,this.project);this.edit(p=>{p.connector={type:'bezier_connector',handleStartM:result.handleStartM,handleEndM:result.handleEndM,curvaturePolicy:'match-neighbors',adjustmentMode:'auto'};});return;}
 }
 private async doImport(input:HTMLInputElement){const file=input.files?.[0];if(!file)return;
 try{if(file.size>2_000_000)throw Error('ファイルサイズ上限は2MBです。');const p=parseProject(await file.text());
     if(!window.confirm('JSONの内容で現在のプロジェクトを置き換えますか？'))return;
     this.editRevision++;this.history.push(structuredClone(this.project));this.project=p;this.plan=buildPlan(p);this.future=[];this.selected=null;this.status='JSON読込完了';this.render();this.map.fit(this.project,this.plan);this.scheduleSave();
   }catch(e){this.notice='JSON読込エラー: '+(e instanceof Error?e.message:String(e));this.renderResults();}
   finally{input.value='';}
 }
}
for(const root of document.querySelectorAll<HTMLElement>('.glp-app-root')){try{new App(root);}catch(e){root.textContent='Glider Landing Planner 起動失敗: '+String(e);console.error(e);}}
