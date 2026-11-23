import test from 'node:test';
import assert from 'node:assert/strict';
import {chartSvg} from '../../assets/js/ui/chart.js';
import {createDemo,trafficPreset} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
test('v0.1.2 chart always renders two independent altitude paths, including disconnected planning',()=>{
 const p=createDemo();p.arrivalLegs=trafficPreset('right');p.connector=null;
 const svg=chartSvg(buildPlan(p),p);
 assert.match(svg,/class="glp-departure-altitude"/);
 assert.match(svg,/class="glp-arrival-altitude"/);
 assert.match(svg,/接続区間（高度線なし）/);
 assert.ok(!svg.includes('class="glp-forwardAll"'));
 assert.match(svg,/高度\(ft, MSL\)/);
});
