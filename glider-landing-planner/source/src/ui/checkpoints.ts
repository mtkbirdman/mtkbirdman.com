import type {Checkpoint,PlanResult,ProjectV1} from '../model/types.js';
import {fmtAlt} from './chart.js';

export interface HeightRow {
 kind:'group'|'height'|'margin';label:string;group?:'departure'|'arrival'|'final';valueM?:number;
}
const groupNames={departure:'離脱側（予想高度）',arrival:'着陸側（逆算必要高度）',final:'ファイナル（目標高度）'};
const encode=(value:string)=>value.replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]??c));
export function landingRelativeHeight(checkpoint:Checkpoint,project:ProjectV1):number {
 return checkpoint.altitudeM-project.landing.elevationMSLM;
}
export function checkpointRows(plan:PlanResult,project:ProjectV1):HeightRow[]{
 const rows:HeightRow[]=[];
 let lastGroup:Checkpoint['group']|null=null;
 for(const cp of plan.checkpoints){
   if(lastGroup!==cp.group){lastGroup=cp.group;rows.push({kind:'group',label:groupNames[cp.group],group:cp.group});}
   if(cp.checkpointKind==='downwind-start')rows.push({kind:'margin',label:'ダウンウィンド開始の高度余裕',valueM:cp.marginM,group:cp.group});
   rows.push({kind:'height',label:cp.label,valueM:landingRelativeHeight(cp,project),group:cp.group});
 }
 return rows;
}
export function checkpointTable(plan:PlanResult,project:ProjectV1,print=false):string{
 const unit=project.display.altitudeUnit;
 const tableClass=print?'glp-print-table':'glp-table';
 const body=checkpointRows(plan,project).map(r=>{
   if(r.kind==='group')return `<tr class="glp-cp-group"><th colspan="2" scope="rowgroup">${encode(r.label)}</th></tr>`;
   const display=r.valueM===undefined?'評価不可':r.kind==='margin'?`${r.valueM>0?'+':''}${fmtAlt(r.valueM,unit)}`:fmtAlt(r.valueM,unit);
   const style=r.kind==='margin'?' class="glp-cp-margin-row"':'';
   return `<tr${style}><td>${encode(r.label)}</td><td>${display}</td></tr>`;
 }).join('');
 return `<table class="${tableClass}"><thead><tr><th>地点</th><th>着陸地点基準高度（${unit}）</th></tr></thead><tbody>${body}</tbody></table>`;
}
