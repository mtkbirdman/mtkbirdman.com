import test from 'node:test';
import assert from 'node:assert/strict';
import {createSinkModel} from '../../assets/js/flight/sink.js';
import {sinkRate} from '../../assets/js/flight/core.js';
const near=(a,b,tol=1e-8)=>assert.ok(Math.abs(a-b)<=tol,`${a} != ${b}`);
const performance=(E=34,V=25)=>({bestLD:E,bestLDAtMps:V});
test('Sink rate retains the exact v0.1.11 quadratic formula at several speeds and banks',()=>{
  const p=performance(), m=createSinkModel(p), E=34, vg=25;
  for (const v of [15,20,25,30,50]) for(const bank of [0,15,30,60]) {
    const expected=v**3/(2*E*vg**2)+vg**2/(2*E*v*Math.cos(bank*Math.PI/180)**2);
    near(m.rate(v,bank),expected,1e-12);
    near(sinkRate(v,bank,p),expected,1e-12);
  }
});

import {simulateLeg} from '../../assets/js/flight/core.js';
import {makeStraight,createDemo} from '../../assets/js/model/defaults.js';

test('a nonfinite sink result leaves geometry evaluable without fabricating altitude',()=>{
  const project=createDemo();
  const state={e:0,n:0,h:500,headingDeg:0,trackDeg:0,bankDeg:0,speedMps:25,groundMps:25,tS:0,sM:0};
  const context={speedMps:25,wind:project.wind,performance:project.performance,sinkModel:{rate:()=>NaN}};
  const result=simulateLeg(state,makeStraight(300),context,1,'departure');
  assert.equal(result.end.altitudeValid,false);
  near(result.end.n,300);
  assert.ok(result.notices.some(n=>n.code==='SINK_RESULT_INVALID'&&n.level==='error'));
});
