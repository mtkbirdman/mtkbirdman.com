import type {ProjectV1,State,PlanResult,Notice,Sample,Checkpoint,FlightLeg} from '../model/types.js';
import {LocalFrame,headingVector,rad,courseOf,add,mul} from '../geometry/geo.js';
import {airGround,headingForTrack,sinkRate,simulatePath} from '../flight/core.js';
import {evaluateConnector} from './connector.js';

export function buildPlan(project:ProjectV1):PlanResult {
 const frame=new LocalFrame(project.landing.point),lp=frame.toLocal(project.landing.point);
 const dep=frame.toLocal(project.departure);
 const ctx={speedMps:project.defaultSpeedMps,wind:project.wind,performance:project.performance};
 const d0=airGround(project.departure.headingDeg,project.defaultSpeedMps,project.wind);
 const departureStart:State={...dep,h:project.departure.altitudeMSLM,headingDeg:project.departure.headingDeg,trackDeg:d0.track,bankDeg:0,speedMps:project.defaultSpeedMps,groundMps:d0.speed,tS:0,sM:0};
 const departure=simulatePath(departureStart,project.departureLegs,ctx,'departure');
 const final=project.landing.final,heading=headingVector(project.landing.runwayTrackDeg);
 const fstart={e:lp.e-heading.e*final.groundLengthM,n:lp.n-heading.n*final.groundLengthM};
 const lastArrivalLeg=project.arrivalLegs[project.arrivalLegs.length-1];
 const lastSpeed=lastArrivalLeg?.speedOverrideMps??project.defaultSpeedMps;
 const solved=headingForTrack(project.landing.runwayTrackDeg,lastSpeed,project.wind);
 const targetHeading=solved?.headingDeg??project.landing.runwayTrackDeg;
 const gd=airGround(targetHeading,lastSpeed,project.wind);
 const finalStart:State={...fstart,h:project.landing.elevationMSLM+final.groundLengthM*Math.tan(rad(final.pathAngleDeg)),headingDeg:targetHeading,trackDeg:project.landing.runwayTrackDeg,bankDeg:0,speedMps:lastSpeed,groundMps:solved?.groundMps??gd.speed,tS:0,sM:0};
 const arrival=simulatePath(finalStart,project.arrivalLegs,ctx,'arrival',true);
 const notices:Notice[]=[...departure.notices,...arrival.notices];
 if(!solved)notices.push({code:'WIND_TRACK_INFEASIBLE',level:'error',message:'ファイナルの対地進入方向を維持できる速度・風の組み合わせではありません。'});
 const connection=project.connector?evaluateConnector(departure.end,arrival.start,project,project.connector.handleStartM,project.connector.handleEndM):null;
 if(connection)notices.push(...connection.notices);else notices.push({code:'CONNECTOR_MISSING',level:'error',message:'離脱側と着陸側が接続されていません。必要高度・到達高度の比較は行いません。'});
 const forwardAll:Sample[]=departure.samples.map(v=>({...v}));
 let predictedB={...departure.end};
 if(connection){forwardAll.push(...connection.samples.slice(1).map(v=>({...v})));predictedB={...connection.end};}
 const predArrivalStart:State={...arrival.start,h:predictedB.h,tS:predictedB.tS,sM:predictedB.sM};
 const predictedArrival=simulatePath(predArrivalStart,project.arrivalLegs,ctx,'arrival',false);
 notices.push(...predictedArrival.notices.filter(n=>n.level==='error'));
 forwardAll.push(...predictedArrival.samples.slice(1).map(v=>({...v})));
 const points=48;
 const finalTarget:Sample[]=[];
 const finalPredicted:Sample[]=[];
 const speedFinal=project.defaultSpeedMps;
 const finalSol=headingForTrack(project.landing.runwayTrackDeg,speedFinal,project.wind);
 const vg=finalSol?.groundMps??airGround(project.landing.runwayTrackDeg,speedFinal,project.wind).speed;
 const finalLoss=vg>0.001?final.groundLengthM*sinkRate(speedFinal,0,project.performance)/vg:0;
 for(let i=0;i<=points;i++){
   const u=i/points,dist=u*final.groundLengthM;
   const e=fstart.e+heading.e*dist,n=fstart.n+heading.n*dist;
   const offsetS=predictedArrival.end.sM+dist,offsetTime=predictedArrival.end.tS+(vg>0?dist/vg:0);
   finalTarget.push({e,n,h:project.landing.elevationMSLM+(final.groundLengthM-dist)*Math.tan(rad(final.pathAngleDeg)),headingDeg:finalSol?.headingDeg??targetHeading,trackDeg:project.landing.runwayTrackDeg,bankDeg:0,speedMps:speedFinal,groundMps:vg,tS:offsetTime,sM:offsetS,segmentId:'final',group:'final'});
   finalPredicted.push({e,n,h:predictedArrival.end.h-(dist/final.groundLengthM)*finalLoss,headingDeg:finalSol?.headingDeg??targetHeading,trackDeg:project.landing.runwayTrackDeg,bankDeg:0,speedMps:speedFinal,groundMps:vg,tS:offsetTime,sM:offsetS,segmentId:'final',group:'final'});
 }
 forwardAll.push(...finalPredicted.slice(1));
 const requiredAll:Sample[]=[];
 const offsetS=predictedB.sM-arrival.start.sM;
 for(const s of arrival.samples){requiredAll.push({...s,sM:s.sM+offsetS});}
 const arrivalNeedEndS=arrival.end.sM+offsetS;
 for(const s of finalTarget){requiredAll.push({...s,sM:arrivalNeedEndS+(s.sM-predictedArrival.end.sM)});}
 const checkpoints:Checkpoint[]=[];
 for(let i=0;i<arrival.segments.length;i++){
   const req=arrival.segments[i].start;
   const pred=predictedArrival.segments[i]?.start;
   if(!pred)continue;
   const leg=arrival.segments[i].leg;
   const label=leg.checkpointLabel??(leg.type==='turn'?'旋回開始 '+(i+1):'レグ開始 '+(i+1));
   checkpoints.push({label,group:'arrival',e:req.e,n:req.n,predictedM:pred.h,requiredM:req.h,marginM:pred.h-req.h});
 }
 if(connection)checkpoints.unshift({label:'自動接続 終点',group:'connector',e:arrival.start.e,n:arrival.start.n,predictedM:predictedB.h,requiredM:arrival.start.h,marginM:predictedB.h-arrival.start.h});
 checkpoints.push({label:'ファイナル開始',group:'final',e:fstart.e,n:fstart.n,predictedM:predictedArrival.end.h,requiredM:finalStart.h,marginM:predictedArrival.end.h-finalStart.h});
 checkpoints.push({label:'進入終端（自然沈下との比較）',group:'final',e:lp.e,n:lp.n,predictedM:finalPredicted[points].h,requiredM:project.landing.elevationMSLM,marginM:finalPredicted[points].h-project.landing.elevationMSLM});
 if(!connection){checkpoints.length=0;forwardAll.splice(departure.samples.length);}
 for(const c of checkpoints){if(c.marginM < -0.5)notices.push({code:'ALTITUDE_DEFICIT',level:'warning',message:`${c.label}: モデル上の高度不足 ${Math.abs(c.marginM/0.3048).toFixed(0)} ft`});}
 notices.push({code:'FINAL_BRAKE_UNKNOWN',level:'info',message:'ファイナル6°は幾何学的な目標線です。ブレーキ・フレア・実際の進入成立性は未評価。'});
 if(project.performance.model==='fixed-ld')notices.push({code:'TURN_SINK_UNMODELED',level:'info',message:'モデルAは旋回追加沈下を未評価。'});
 for(const l of [...project.departureLegs,...project.arrivalLegs]){if(l.type==='turn'&&l.bankDeg>project.maxBankDeg){notices.push({code:'BANK_LIMIT_EXCEEDED',level:'error',message:`レグ ${l.id} の指定バンク角 ${l.bankDeg}° が上限 ${project.maxBankDeg}° を超えています。`,legId:l.id});}}
 if(connection&&Math.abs(connection.end.headingDeg-arrival.start.headingDeg)>0.3){notices.push({code:'BOUNDARY_DISCONTINUITY',level:'warning',message:'Bezierと着陸側経路の機首方位が不連続になる場合があります。'});}
 const valid=!notices.some(n=>n.level==='error');
 return {departure,arrival,connector:connection,forwardAll,requiredAll,finalTarget,checkpoints,notices,arrivalStart:arrival.start,finalStart,landing:lp,valid};
}
