import type {PlanResult,ProjectV1,Sample} from '../model/types.js';
export const alt=(m:number,u:'m'|'ft')=>u==='ft'?m/0.3048:m;
export const fmtAlt=(m:number,u:'m'|'ft')=>`${Math.round(alt(m,u)).toLocaleString('ja-JP')} ${u}`;
export const fmtElapsed=(seconds:number):string=>{
 const s=Math.max(0,Math.round(seconds));
 return `${Math.floor(s/60)}:${String(s%60).padStart(2,'0')}`;
};

/**
 * Use genuine elapsed seconds, not ENU distance or a distance-based display offset.
 * Reverse-integrated arrival samples carry negative tS; shift them to the actual
 * end of the connector, without modifying the canonical flight/plan data.
 * The two curves intentionally remain separate even when their x coordinates meet.
 */
export function altitudeTimeSeries(plan:PlanResult):{
 departure:Sample[];arrival:Sample[];durationS:number;timingProvisional:boolean;
}{
 const finite=(s:Sample)=>Number.isFinite(s.tS)&&Number.isFinite(s.h);
 const departure=[...plan.departure.samples,...(plan.connector?.samples.slice(1)??[])].filter(finite);
 const joinedTime=departure.at(-1)?.tS??0;
 const arrivalStartRaw=plan.arrival.samples[0]?.tS??plan.arrival.start.tS;
 const arrivalShift=joinedTime-arrivalStartRaw;
 const arrival=plan.arrival.samples.filter(finite).map(s=>({...s,tS:s.tS+arrivalShift}));
 const arrivalEnd=arrival.at(-1)?.tS??joinedTime;
 const finalStartRaw=plan.finalTarget[0]?.tS??0;
 arrival.push(...plan.finalTarget.slice(1).filter(finite).map(s=>({...s,tS:arrivalEnd+s.tS-finalStartRaw})));
 let durationS=0;for(const s of departure)durationS=Math.max(durationS,s.tS);for(const s of arrival)durationS=Math.max(durationS,s.tS);
 return {departure,arrival,durationS,timingProvisional:!plan.connector||!plan.connector.feasible};
}

export function chartSvg(plan:PlanResult,project:ProjectV1,width=780,height=200):string{
 // Blue: forecast including quintic Bezier connector. Green: independent arrival requirement.
 const {departure,arrival,durationS,timingProvisional}=altitudeTimeSeries(plan);
 const all=[...departure,...arrival];if(!all.length)return '<svg viewBox="0 0 780 200"></svg>';
 const left=52,right=16,top=23,bottom=35;
 let ymin=Infinity,ymax=-Infinity;for(const v of all){if(v.h<ymin)ymin=v.h;if(v.h>ymax)ymax=v.h;}
 const ypad=Math.max(20,(ymax-ymin)*.14),y0=ymin-ypad,y1=ymax+ypad;
 const px=(t:number)=>left+Math.max(0,t)/Math.max(1,durationS)*(width-left-right);
 const py=(y:number)=>top+(y1-y)/Math.max(1,y1-y0)*(height-top-bottom);
 const trace=(pts:Sample[],color:string,cssClass:string)=>{
   if(!pts.length)return '';
   const stride=Math.max(1,Math.ceil(pts.length/1000));
   const subset=pts.filter((_,i)=>i%stride===0||i===pts.length-1);
   const path=subset.map((s,i)=>`${i?'L':'M'}${px(s.tS).toFixed(1)},${py(s.h).toFixed(1)}`).join(' ');
   return `<path class="${cssClass}" d="${path}" fill="none" stroke="${color}" stroke-width="2.3" stroke-linejoin="round"/>`;
 };
 const u=project.display.altitudeUnit;
 let ticks='';for(let i=0;i<=4;i++){
   const y=y0+(y1-y0)*i/4,yp=py(y);
   ticks+=`<line x1="${left}" y1="${yp.toFixed(1)}" x2="${width-right}" y2="${yp.toFixed(1)}" stroke="#dbe4ef" stroke-width=".8"/>`;
   ticks+=`<text x="${left-8}" y="${(yp+4).toFixed(1)}" text-anchor="end" font-size="11" fill="#607086">${Math.round(alt(y,u))}</text>`;
 }
 // Horizontal axis is elapsed flight time, with actual time labels and no distance offset.
 const stepCount=4;
 for(let i=0;i<=stepCount;i++){
   const t=durationS*i/stepCount,x=px(t);
   ticks+=`<line x1="${x.toFixed(1)}" y1="${height-bottom}" x2="${x.toFixed(1)}" y2="${height-bottom+4}" stroke="#94a3b8"/>`;
   ticks+=`<text x="${x.toFixed(1)}" y="${height-17}" font-size="10" text-anchor="middle" fill="#607086">${fmtElapsed(t)}</text>`;
 }
 const label=timingProvisional?'未接続・制約未成立：着陸側時刻は仮配置':'離脱からの推定経過飛行時間';
 return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="離脱側予想高度（Bezier接続を含む）と着陸側逆算必要高度。横軸は飛行時間。">
 <rect x="0" y="0" width="${width}" height="${height}" fill="#fff"/><text x="9" y="12" font-size="11" fill="#64748b">高度(${u}, MSL)</text>${ticks}
 ${trace(departure,'#2563eb','glp-departure-altitude')}${trace(arrival,'#059669','glp-arrival-altitude')}
 <text x="${left}" y="${height-2}" font-size="10" fill="#64748b">${label}</text><text x="${width-right}" y="${height-2}" font-size="10" fill="#64748b" text-anchor="end">計画時間 ${fmtElapsed(durationS)}</text>
 </svg>`;
}
export function printMapSvg(plan:PlanResult,width=780,height=420):string {
 const groups:{points:Sample[];color:string;dashed?:boolean}[]=[
  {points:plan.departure.samples,color:'#2563eb'},
  {points:plan.connector?.samples??[],color:'#8b5cf6'},
  {points:plan.arrival.samples,color:'#059669'},
  {points:plan.finalTarget,color:'#059669',dashed:true}
 ];
 const all=groups.flatMap(g=>g.points).filter(s=>Number.isFinite(s.e)&&Number.isFinite(s.n));if(all.length===0)return '';
 const xmin=Math.min(...all.map(p=>p.e)),xmax=Math.max(...all.map(p=>p.e)),ymin=Math.min(...all.map(p=>p.n)),ymax=Math.max(...all.map(p=>p.n));
 const margin=40,scale=Math.min((width-margin*2)/Math.max(50,xmax-xmin),(height-margin*2)/Math.max(50,ymax-ymin));
 const cx=(xmin+xmax)/2,cy=(ymin+ymax)/2;
 const px=(e:number)=>width/2+(e-cx)*scale,py=(n:number)=>height/2-(n-cy)*scale;
 const traces=groups.map(g=>{
  if(g.points.length<2)return '';
  const stride=Math.max(1,Math.ceil(g.points.length/600));
  const points=g.points.filter((_,i)=>i%stride===0||i===g.points.length-1);
  const path=points.map((s,i)=>`${i?'L':'M'}${px(s.e).toFixed(1)},${py(s.n).toFixed(1)}`).join(' ');
  return `<path d="${path}" fill="none" stroke="${g.color}" stroke-width="3" ${g.dashed?'stroke-dasharray="6 5"':''}/>`;
 }).join('');
 const d=plan.departure.start,l=plan.landing;
 return `<svg viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" role="img" aria-label="経路図（背景地図なし）">
 <rect width="${width}" height="${height}" fill="#f6f9fc"/><path d="M ${width/2} 0 V ${height} M 0 ${height/2} H ${width}" stroke="#cbd5e1" stroke-dasharray="3 6"/>
 ${traces}
 <circle cx="${px(d.e)}" cy="${py(d.n)}" r="5" fill="#2563eb"/><circle cx="${px(l.e)}" cy="${py(l.n)}" r="5" fill="#059669"/>
 <text x="${px(d.e)+7}" y="${py(d.n)-8}" font-size="12" fill="#1e293b">離脱</text><text x="${px(l.e)+7}" y="${py(l.n)-8}" font-size="12" fill="#1e293b">着陸</text>
 <path d="M ${width-35} 42 V 16 M ${width-39} 22 L ${width-35} 16 L ${width-31} 22" stroke="#334155" stroke-width="2" fill="none"/><text x="${width-39}" y="12" font-size="12">N</text>
 <text x="12" y="${height-10}" fill="#64748b" font-size="11">ENU座標による経路図（背景地図なし・縮尺は図枠に合わせて調整）</text></svg>`;
}
