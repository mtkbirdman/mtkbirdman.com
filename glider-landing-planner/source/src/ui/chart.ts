import type {PlanResult,ProjectV1,Sample} from '../model/types.js';
export const alt=(m:number,u:'m'|'ft')=>u==='ft'?m/0.3048:m;
export const fmtAlt=(m:number,u:'m'|'ft')=>`${Math.round(alt(m,u)).toLocaleString('ja-JP')} ${u}`;
export function chartSvg(plan:PlanResult,project:ProjectV1,width=780,height=200):string{
 // Render two INDEPENDENT altitude functions: departure forecast vs arrival required height.
 // Deliberately do not draw a predicted connecting profile or bridge the Bezier segment.
 const departure=plan.departure.samples.filter(s=>Number.isFinite(s.sM)&&Number.isFinite(s.h));
 const reqStartX=plan.departure.end.sM+(plan.connector?.lengthM??Math.max(50,Math.hypot(plan.arrival.start.e-plan.departure.end.e,plan.arrival.start.n-plan.departure.end.n)));
 const shift=reqStartX-plan.arrival.start.sM;
 const arrival=[...plan.arrival.samples,...plan.finalTarget.slice(1)].filter(s=>Number.isFinite(s.sM)&&Number.isFinite(s.h)).map(s=>({...s,sM:s.sM+shift}));
 const all=[...departure,...arrival];if(!all.length)return '<svg viewBox="0 0 780 200"></svg>';
 const left=48,right=14,top=15,bottom=28;
 const xmin=Math.min(...all.map(v=>v.sM)),xmax=Math.max(...all.map(v=>v.sM)),ymin=Math.min(...all.map(v=>v.h)),ymax=Math.max(...all.map(v=>v.h));
 const ypad=Math.max(20,(ymax-ymin)*.14),y0=ymin-ypad,y1=ymax+ypad;
 const px=(x:number)=>left+(x-xmin)/Math.max(1,xmax-xmin)*(width-left-right);
 const py=(y:number)=>top+(y1-y)/Math.max(1,y1-y0)*(height-top-bottom);
 const trace=(pts:Sample[],color:string,cssClass:string)=>{
   if(!pts.length)return '';
   const stride=Math.max(1,Math.ceil(pts.length/700));
   const subset=pts.filter((_,i)=>i%stride===0||i===pts.length-1);
   const path=subset.map((s,i)=>`${i?'L':'M'}${px(s.sM).toFixed(1)},${py(s.h).toFixed(1)}`).join(' ');
   return `<path class="${cssClass}" d="${path}" fill="none" stroke="${color}" stroke-width="2.3" stroke-linejoin="round"/>`;
 };
 const u=project.display.altitudeUnit;
 let ticks='';for(let i=0;i<=4;i++){
   const y=y0+(y1-y0)*i/4,yp=py(y);
   ticks+=`<line x1="${left}" y1="${yp.toFixed(1)}" x2="${width-right}" y2="${yp.toFixed(1)}" stroke="#dbe4ef" stroke-width=".8"/>`;
   ticks+=`<text x="${left-8}" y="${(yp+4).toFixed(1)}" text-anchor="end" font-size="11" fill="#607086">${Math.round(alt(y,u))}</text>`;
 }
 const gapBegin=px(plan.departure.end.sM),gapEnd=px(reqStartX);
 const gapBand=gapEnd-gapBegin>8?`<rect class="glp-chart-break" x="${gapBegin.toFixed(1)}" y="${top}" width="${(gapEnd-gapBegin).toFixed(1)}" height="${height-bottom-top}" fill="#94a3b8" fill-opacity=".075"/><text x="${((gapBegin+gapEnd)/2).toFixed(1)}" y="${top+13}" text-anchor="middle" font-size="10" fill="#64748b">接続区間（高度線なし）</text>`:'';
 return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="離脱側予想高度と着陸側逆算必要高度（不連続）">
 <rect x="0" y="0" width="${width}" height="${height}" fill="#fff"/><text x="9" y="12" font-size="11" fill="#64748b">高度(${u}, MSL)</text>${ticks}${gapBand}
 ${trace(departure,'#2563eb','glp-departure-altitude')}${trace(arrival,'#059669','glp-arrival-altitude')}
 <text x="${left}" y="${height-6}" font-size="11" fill="#64748b">離脱</text><text x="${width-right}" y="${height-6}" font-size="11" fill="#64748b" text-anchor="end">累積計画経路距離 ${Math.round(xmax)} m</text>
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
