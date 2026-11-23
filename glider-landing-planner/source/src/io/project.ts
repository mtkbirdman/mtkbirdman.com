import type {ProjectV1,FlightLeg} from '../model/types.js';
const between=(n:unknown,min:number,max:number)=>typeof n==='number'&&Number.isFinite(n)&&n>=min&&n<=max;
const obj=(v:unknown):Record<string,unknown>=>{if(!v||typeof v!=='object'||Array.isArray(v))throw Error('オブジェクトが必要です。');return v as Record<string,unknown>};
function check(ok:boolean,msg:string){if(!ok)throw Error(msg);}
function validateLegs(legs:unknown,ids:Set<string>,label:string){check(Array.isArray(legs)&&legs.length<=100,label+' は100レグ以下で指定してください。');
 for(const input of legs as unknown[]){const l=obj(input),id=l.id;check(typeof id==='string'&&id.length>0&&id.length<=80&&!ids.has(id),label+': レグIDが不正または重複しています。');ids.add(id as string);
 check(['heading','ground-track'].includes(String(l.guidance)),label+': 風補正モードが不正です。');
 check(l.speedOverrideMps===null||between(l.speedOverrideMps,5,100),label+': 速度上書きが不正です。');
 check(l.directionDeg===null||l.directionDeg===undefined||between(l.directionDeg,0,360),label+': 方位指定が不正です。');
 if(l.checkpointLabel!==undefined)check(typeof l.checkpointLabel==='string'&&l.checkpointLabel.length<=70,label+': チェックポイント名が長すぎます。');
 if(l.checkpointKind!==undefined)check(l.checkpointKind==='downwind-start',label+': チェックポイント種別が不正です。');
 if(l.type==='straight')check(between(l.lengthM,1,50000),label+': 直線距離が不正です。');
 else if(l.type==='turn'){check(['left','right'].includes(String(l.direction)),label+': 旋回方向が不正です。');check(between(l.angleDeg,1,720)&&between(l.bankDeg,1,70),label+': 旋回角・バンク角が不正です。');}
 else throw Error(label+': 未対応のレグ種別です。');
 }
}
export function validateProject(input:unknown):ProjectV1 {
 const p=obj(input);
 check(p.schemaVersion==='1.0','未対応のschemaVersionです（対応: 1.0）。');
 check(typeof p.engineVersion==='string'&&p.engineVersion.length<70,'engineVersionが不正です。');
 check(typeof p.projectId==='string'&&p.projectId.length<=100,'projectIdが不正です。');
 check(typeof p.name==='string'&&p.name.length<=120,'計画名は120文字以内です。');
 check(between(p.defaultSpeedMps,5,100),'基本速度は5〜100 m/sの範囲にしてください。');
 check(between(p.maxBankDeg,5,70),'最大バンク角は5〜70°の範囲にしてください。');
 const perf=obj(p.performance);check(['fixed-ld','quadratic-drag'].includes(String(perf.model)),'計算モデルが不正です。');
 check(between(perf.fixedLD,2,100)&&between(perf.bestLD,2,100)&&between(perf.bestLDAtMps,5,100),'性能パラメータが不正です。');
 const wind=obj(p.wind);check(between(wind.fromDeg,0,360)&&between(wind.speedMps,0,80),'風向・風速が不正です。');
 const d=obj(p.departure);check(between(d.latDeg,-90,90)&&between(d.lonDeg,-180,180),'離脱地点の緯度経度が不正です。');
 check(between(d.altitudeMSLM,-500,15000)&&between(d.headingDeg,0,360),'離脱高度・機首方位が不正です。');
 const l=obj(p.landing),lp=obj(l.point),final=obj(l.final);
 check(between(lp.latDeg,-90,90)&&between(lp.lonDeg,-180,180),'着陸地点の座標が不正です。');
 check(between(l.elevationMSLM,-500,15000)&&between(l.runwayTrackDeg,0,360),'着陸標高・進入方位が不正です。');
 check(between(final.pathAngleDeg,0.1,25)&&between(final.groundLengthM,20,10000),'ファイナル角度・距離が不正です。');
 const ids=new Set<string>();validateLegs(p.departureLegs,ids,'離脱側レグ');validateLegs(p.arrivalLegs,ids,'着陸側レグ');
 // v0.1.1 compatibility: migrate the legacy Japanese label once, retaining its identity when renamed later.
 const arrival=p.arrivalLegs as FlightLeg[];if(!arrival.some(v=>v.checkpointKind==='downwind-start')){
   const legacy=arrival.find(v=>v.checkpointLabel==='ダウンウィンド開始');if(legacy)legacy.checkpointKind='downwind-start';
 }
 if(p.connector!==null){const c=obj(p.connector);check(c.type==='bezier_connector'&&c.curvaturePolicy==='match-neighbors','接続曲線種別が不正です。');check(between(c.handleStartM,1,50000)&&between(c.handleEndM,1,50000),'Bezierハンドル長が不正です。');}
 const disp=obj(p.display);check(['m','ft'].includes(String(disp.altitudeUnit)),'高度単位が不正です。');
 for(const k of ['showDistanceRings','showBearingLines','showWindArrow','showHandles'])check(typeof disp[k]==='boolean','表示設定 '+k+' が不正です。');
 // Additive migration; preserve schemaVersion 1.0 to read existing saved files and IndexedDB snapshots.
 const defaults:Record<string,boolean|number>={showDownwindMargin:true,showCheckpointAltitudes:false,distanceRingOpacity:0.5,bearingLineOpacity:0.5};
 for(const [k,v] of Object.entries(defaults))if(disp[k]===undefined)disp[k]=v;
 for(const k of ['showDownwindMargin','showCheckpointAltitudes'])check(typeof disp[k]==='boolean','表示設定 '+k+' が不正です。');
 for(const k of ['distanceRingOpacity','bearingLineOpacity'])check(between(disp[k],0,1),'不透明度 '+k+' は0〜1で指定してください。');
 return input as ProjectV1;
}
export function parseProject(text:string):ProjectV1{if(text.length>2_000_000)throw Error('ファイルが大きすぎます（上限2MB）。');
 const candidate=JSON.parse(text) as unknown;return validateProject(candidate);}
export function serializeProject(project:ProjectV1):string {return JSON.stringify(validateProject(project),null,2);}
export function downloadProject(project:ProjectV1){const blob=new Blob([serializeProject(project)],{type:'application/json'});const url=URL.createObjectURL(blob);const a=document.createElement('a');a.href=url;a.download=(project.name.replace(/[^\p{L}\p{N}_-]+/gu,'_')||'glider_plan')+'.glp.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
