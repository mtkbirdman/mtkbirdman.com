import test from 'node:test';
import assert from 'node:assert/strict';
import {chartSvg,altitudeTimeSeries,fmtElapsed} from '../../assets/js/ui/chart.js';
import {createDemo,trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';

const near=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<tol,`${a} != ${b}`);
test('v0.1.3 altitude chart has independent blue/green paths and time ticks',()=>{
 const p=createDemo();p.arrivalLegs=trafficPreset('right');p.connector=null;
 const svg=chartSvg(buildPlan(p),p);
 assert.match(svg,/class="glp-departure-altitude"/);
 assert.match(svg,/class="glp-arrival-altitude"/);
 assert.doesNotMatch(svg,/接続区間（高度線なし）/);
 assert.doesNotMatch(svg,/累積計画経路距離/);
 assert.match(svg,/未接続・制約未成立：着陸側時刻は仮配置/);
 assert.match(svg,/計画時間 \d+:\d\d/);
 assert.match(svg,/高度\(ft, MSL\)/);
 assert.match(svg,/>0:00</);
});
test('v0.1.3 connector time and altitude are included in departure blue line, no cumulative double offset',()=>{
 const p=createDemo();p.arrivalLegs=trafficPreset('left');p.wind.speedMps=5;
 const plan=buildPlan(p), series=altitudeTimeSeries(plan);
 assert.ok(plan.connector);
 assert.equal(series.departure.length,plan.departure.samples.length+plan.connector.samples.length-1);
 near(series.departure.at(-1).tS,plan.connector.end.tS);
 near(series.departure.at(-1).h,plan.connector.end.h);
 near(series.arrival[0].tS,series.departure.at(-1).tS);
 const reverseDuration=plan.arrival.end.tS-plan.arrival.start.tS;
 const finalDuration=plan.finalTarget.at(-1).tS-plan.finalTarget[0].tS;
 near(series.durationS,plan.connector.end.tS+reverseDuration+finalDuration);
 near(series.arrival.at(-1).tS,series.durationS);
 assert.ok(series.arrival.every((s,i)=>!i||s.tS>=series.arrival[i-1].tS-1e-8));
 assert.equal(series.arrival.at(-1).h,p.landing.elevationMSLM);
 assert.match(chartSvg(plan,p),/class="glp-departure-altitude"/);
});
test('v0.1.3 time formatting rounds seconds and displays min:sec',()=>{
 assert.equal(fmtElapsed(0),'0:00');assert.equal(fmtElapsed(61.8),'1:02');
 assert.equal(fmtElapsed(3600),'60:00');
});
