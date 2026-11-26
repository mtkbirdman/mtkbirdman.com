import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo,makeStraight,makeTurn,trafficPreset,ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {parseProject,serializeProject} from '../../assets/js/io/project.js';
const read=(file)=>readFileSync(new URL(file,import.meta.url),'utf8');

test('v0.1.8 WordPress version and engine version match',()=>{
 assert.equal(ENGINE_VERSION,'0.1.20-test');
 assert.match(read('../../glider-landing-planner.php'),/Version: 0\.1\.20/);
 assert.match(read('../src/main.ts'),/テスト版 v0\.1\.20/);
});
test('v0.1.7 departure checkpoints count turns and straights independently',()=>{
 const p=createDemo();p.departureLegs=[makeStraight(200),makeTurn('left',90),makeStraight(150),makeTurn('right',270)];
 const labels=buildPlan(p).checkpoints.filter(c=>c.group==='departure').map(c=>c.label);
 assert.deepEqual(labels,['離脱地点','離脱側 直線1 終点','離脱側 第1旋回 終点','離脱側 直線2 終点','離脱側 第2旋回 終点']);
});
test('v0.1.7 arrival checkpoints count only unnamed added legs; preset labels unchanged',()=>{
 const p=createDemo();p.arrivalLegs=[makeTurn('left',120),makeStraight(250),...trafficPreset('left'),makeTurn('right',180),makeStraight(350)];
 const labels=buildPlan(p).checkpoints.filter(c=>c.kind==='arrival-start').map(c=>c.label);
 assert.deepEqual(labels,['着陸側 第1旋回 開始','着陸側 直線1 開始','ダウンウィンド開始','ベースターン開始','ベースレグ開始','ファイナルターン開始','着陸側 第2旋回 開始','着陸側 直線2 開始']);
});
test('v0.1.7 user-customized checkpoint names survive unchanged',()=>{
 const p=previouslyPlannedProject();p.departureLegs[0].checkpointLabel='訓練ポイントA';p.arrivalLegs=[{...makeTurn(),checkpointLabel:'手動ポイントB'}];
 const names=buildPlan(p).checkpoints.map(c=>c.label);
 assert.ok(names.includes('訓練ポイントA'));assert.ok(names.includes('手動ポイントB'));
});
test('v0.1.7 main reuses one editor for sidebar and bottom-sheet',()=>{
 const main=read('../src/main.ts');
 assert.match(main,/private legEditor\(side:\s*Side,\s*leg:\s*FlightLeg,\s*i:\s*number\)/);
 assert.match(main,/this\.legEditor\(side,\s*leg,\s*i\)/);
 assert.match(main,/this\.legEditor\(side,\s*leg,\s*index\)/);
 assert.match(main,/leg\.type\s*===\s*'turn'\s*\?\s*'旋回速度'\s*:\s*'滑空速度'/);
 assert.match(main,/onSelectLeg:\s*\(id\)\s*=>\s*this\.selectLeg\(id,\s*true\)/);
 assert.match(main,/onSelectPoint:\s*\(point\)\s*=>\s*\{\s*this\.selected\s*=\s*null/);
 assert.match(main,/this\.openMobilePanel\(\{\s*kind:\s*'point',\s*point\s*\}\)/);
 assert.match(main,/onAddAtEnd:\s*\(side, type\)\s*=>/);
 assert.match(main,/this\.openMobilePanel\(\{\s*kind:\s*'add',\s*side\s*\}\)/);
 assert.match(main,/data-action="preset-left"/);
 assert.match(main,/this\.selectLeg\(leg\.id,\s*true\)/);
});
test('v0.1.7 5th-order Bezier connector policy remains unchanged',()=>{
 const p=previouslyPlannedProject();assert.equal(p.connector.curvaturePolicy,'match-neighbors');
 const src=read('../src/planning/connector.ts');
 assert.match(src,/curvature/);
 assert.doesNotMatch(read('../src/main.ts'),/rollout|ロールアウト処理/);
});
test('v0.1.7 mobile style and print hide editing controls',()=>{
 const style=read('../../assets/css/app.css'),print=read('../../assets/css/print.css');
 assert.match(style,/\.glp-mobile-editor:not\(\[hidden\]\)/);
 assert.match(style,/orientation:portrait/);
 assert.match(style,/\.glp-map-add-badge/);
 assert.match(print,/\.glp-map-add-point/);
});
test('v0.1.7 old project data and checkpoint labels remain readable',()=>{
 const p=createDemo();p.engineVersion='0.1.6-test';p.departureLegs=[makeTurn('right',1080),makeStraight(100)];p.arrivalLegs=trafficPreset('right');
 const data=parseProject(serializeProject(p));assert.equal(data.departureLegs[0].angleDeg,1080);assert.equal(data.arrivalLegs[0].checkpointLabel,'ダウンウィンド開始');
});
