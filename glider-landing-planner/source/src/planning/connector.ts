import type {State,Point2,ConnectionResult,Notice,Sample,ProjectV1} from '../model/types.js';
import {headingVector,courseOf,rad,deg,dot,magnitude,mul,add,sub,leftNormal} from '../geometry/geo.js';
import {point,derivative,curvature,controlPoints,selfIntersects} from '../geometry/bezier.js';
import {windVector,headingForTrack,sinkRate,GRAVITY} from '../flight/core.js';
const EPS=1e-8;
export function endpointCurvature(state:State,bankDeg:number,project:ProjectV1):number {
 if(Math.abs(bankDeg)<0.01)return 0;
 const v=state.speedMps,omega=GRAVITY*Math.tan(rad(bankDeg))/v;
 const Vg=add(mul(headingVector(state.headingDeg),v),windVector(project.wind));
 const dV={e:v*Math.cos(rad(state.headingDeg))*omega,n:-v*Math.sin(rad(state.headingDeg))*omega};
 return (Vg.e*dV.n-Vg.n*dV.e)/(magnitude(Vg)**3);
}
interface ConnectorGeometry {points:Point2[];cp:Point2[];lengthM:number;maxBankDeg:number;invalid:boolean;notices:Notice[];samples:Sample[];end:State;lossM:number;}
export function evaluateConnector(start:State,endBoundary:State,project:ProjectV1,a:number,b:number,skipHeavy=false):ConnectionResult {
 const A:Point2={e:start.e,n:start.n},B:Point2={e:endBoundary.e,n:endBoundary.n};
 const tA=headingVector(start.trackDeg),tB=headingVector(endBoundary.trackDeg);
 const kA=endpointCurvature(start,start.bankDeg,project),kB=endpointCurvature(endBoundary,endBoundary.bankDeg,project);
 const cp=controlPoints(A,B,tA,tB,kA,kB,a,b);
 const approximateSpan=magnitude(sub(B,A))+a+b;
 const steps=skipHeavy?56:Math.min(2000,Math.max(160,Math.ceil(approximateSpan/3)));
 const points:Point2[]=[],samples:Sample[]=[],notices:Notice[]=[];
 let invalid=false,maxBank=0,length=0,loss=0,time=0;
 const w=windVector(project.wind),V=start.speedMps;
 let prev:{p:Point2;g:number;sink:number}|null=null;
 let nowState={...start};
 let prevTrack=start.trackDeg;
 for(let i=0;i<=steps;i++){
   const u=i/steps,p=point(cp,u),d=derivative(cp,u),dsdu=magnitude(d),t=dsdu>EPS?mul(d,1/dsdu):tA;
   const k=dsdu>EPS?curvature(cp,u):NaN;
   let G=0,heading=prevTrack,bank=0;
   const Wt=dot(w,t),Wn=dot(w,leftNormal(t)),q2=V*V-Wn*Wn;
   if(!Number.isFinite(k)||dsdu<1e-5||q2<=EPS){invalid=true;}
   else{
     const q=Math.sqrt(q2);G=Wt+q;
     if(G<=EPS){invalid=true;}else{
       const air=sub(mul(t,G),w);heading=courseOf(air.e,air.n);
       const omega=-k*G*G/q;
       bank=deg(Math.atan(omega*V/GRAVITY));
       maxBank=Math.max(maxBank,Math.abs(bank));
     }
   }
   points.push(p);
   const ssink=sinkRate(V,bank,project.performance);
   if(prev){const ds=magnitude(sub(p,prev.p));length+=ds;
     if(G>EPS&&prev.g>EPS){const dt=ds*0.5*(1/G+1/prev.g);time+=dt;const dh=dt*0.5*(ssink+prev.sink);loss+=dh;}
     else invalid=true;
   }
   const s:Sample={e:p.e,n:p.n,h:start.h-loss,headingDeg:heading,trackDeg:courseOf(t.e,t.n),bankDeg:bank,speedMps:V,groundMps:G,tS:start.tS+time,sM:start.sM+length,segmentId:'connector',group:'connector'};
   samples.push(s);prev={p,g:G,sink:ssink};prevTrack=s.trackDeg;
 }
 if(selfIntersects(points)){invalid=true;notices.push({code:'CONNECTOR_SELF_INTERSECTION',level:'error',message:'接続曲線が自己交差しています。'});}
 if(maxBank>project.maxBankDeg+1e-6){invalid=true;notices.push({code:'BANK_LIMIT_EXCEEDED',level:'error',message:`Bezierの必要バンク角（最大${maxBank.toFixed(1)}°）が設定上限を超過。`});}
 if(invalid&&!notices.some(n=>n.code==='CONNECTOR_SELF_INTERSECTION'))notices.push({code:'CONNECTOR_NO_SOLUTION',level:'error',message:'曲線追従に必要な速度・曲率条件を満たしません。'});
 const end= samples.length ? {...samples[samples.length-1]}: {...start};
 if(samples.length){
   const startHeading=samples[0].headingDeg,endHeading=samples[samples.length-1].headingDeg;
   const diffA=Math.abs(((startHeading-start.headingDeg+540)%360)-180);
   const diffB=Math.abs(((endHeading-endBoundary.headingDeg+540)%360)-180);
   if(diffA>0.3||diffB>0.3)notices.push({code:'BOUNDARY_DISCONTINUITY',level:'warning',message:`接続部の機首方位差: 始端 ${diffA.toFixed(1)}° / 終端 ${diffB.toFixed(1)}°。`});
   if(Math.abs(endBoundary.speedMps-V)>1e-4)notices.push({code:'SPEED_STEP_IDEALIZED',level:'warning',message:'Bezier終点で個別速度との不連続があります。'});
 }
 if(project.performance.model==='fixed-ld')notices.push({code:'TURN_SINK_UNMODELED',level:'info',message:'モデルAのBezier接続区間では旋回追加沈下を評価していません。'});
 return {points,samples,end,lengthM:length,lossM:loss,maxBankDeg:maxBank,feasible:!invalid,notices,handleStartM:a,handleEndM:b};
}
export function chooseConnector(start:State,end:State,project:ProjectV1):ConnectionResult {
 const length=Math.hypot(end.e-start.e,end.n-start.n);
 if(length<1){return evaluateConnector(start,end,project,5,5);}
 let best:ConnectionResult|null=null,bestScore=Infinity;
 for(const x of [0.16,0.27,0.40,0.58,0.8,1.1])for(const y of [0.16,0.27,0.40,0.58,0.8,1.1]){
   const r=evaluateConnector(start,end,project,Math.max(10,length*x),Math.max(10,length*y),true);
   const score=(r.feasible?0:1e7)+r.maxBankDeg*120+r.lengthM+Math.max(0,r.maxBankDeg-project.maxBankDeg)*1e4;
   if(score<bestScore){bestScore=score;best=r;}
 }
 if(!best)throw Error('connector search failed');
 const r=evaluateConnector(start,end,project,best.handleStartM,best.handleEndM);
 if(!r.feasible)r.notices.push({code:'CONNECTOR_NO_SOLUTION',level:'error',message:'探索範囲内で制約を満たすBezier接続は見つかりませんでした。'});
 return r;
}
