import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync} from 'node:fs';
import {createDemo,ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {checkpointTable} from '../../assets/js/ui/checkpoints.js';
const read=(p)=>readFileSync(new URL(p,import.meta.url),'utf8');

test('v0.1.5 default limit is 60 degrees; existing record remains user-editable',()=>{
 const p=createDemo();assert.equal(ENGINE_VERSION,'0.1.14-test');assert.equal(p.maxBankDeg,60);
 p.maxBankDeg=45;assert.equal(buildPlan(p).notices.find(n=>n.code==='invalid-maximum-bank')?.level,undefined);
});

test('v0.1.5 checkpoint rows carry stable leg ids for cross-highlighting',()=>{
 const p=createDemo();p.arrivalLegs.push({id:'chosen-123',type:'straight',lengthM:200,guidance:'heading',directionDeg:null,speedOverrideMps:null});
 const html=checkpointTable(buildPlan(p),p);
 assert.match(html,/data-leg-id="chosen-123"/);
 const main=read('../src/main.ts');
 assert.match(main,/updateCheckpointSelection\(\)/);
 assert.match(main,/row\.dataset\.legId\s*===\s*this\.selected/);
});

test('v0.1.7 arrival add controls precede arrival list, final length is at bottom',()=>{
 const main=read('../src/main.ts');
 assert.match(main,/side\s*===\s*'arrival'\s*\?\s*adds\s*\+\s*list\s*:\s*list\s*\+\s*adds/);
 assert.match(main,/右場周プリセット<\/button><\/div>\$\{this\.legList\('arrival'\)/);
 assert.match(main,/新しいレグは最上流へ追加され、一覧は実際の飛行順です。<\/p>\$\{this\.control\('landing\.final\.groundLengthM'/);
 assert.doesNotMatch(main,/風向（吹いてくる）/);
 assert.match(main,/const preservedScroll\s*=\s*sidebarScroll\.scrollTop/);
 assert.match(main,/sidebarScroll\.scrollTop\s*=\s*preservedScroll/);
});

test('v0.1.5 print preview clones live Leaflet map, preserves GSI tiles and all sections',()=>{
 const print=read('../src/ui/print.ts');
 const css=read('../../assets/css/print.css');
 assert.match(print,/mapEl\.cloneNode\(true\)/);
 assert.match(print,/leaflet-tile-pane img\.leaflet-tile/);
 assert.match(print,/glp-report-map-frame/);
 for(const title of ['飛行経路図','高度プロファイル','計算結果','飛行条件','チェックポイント高度','計算上の診断'])assert.ok(print.includes(title));
 assert.match(print,/window\.open\(\s*''\s*,\s*'_blank'\s*\)/);
 assert.match(print,/popup\.print\(\)/);
 assert.doesNotMatch(css,/body \*\s*\{\s*visibility:hidden/);
 assert.doesNotMatch(css,/float:\s*(left|right)/);
 assert.match(css,/@page\{size:A4 landscape/);
 assert.match(css,/\.glp-report-chart/);
});

test('v0.1.5 popup stylesheet is independent of Cocoon page print layout',()=>{
 const print=read('../src/ui/print.ts');
 const css=read('../../assets/css/print.css');
 assert.match(print,/new URL\('\.\.\/\.\.\/css\/print\.css',\s*import\.meta\.url\)/);
 assert.match(print,/国土地理院/);
 assert.match(css,/\.glp-report-document\{width:/);
 assert.match(css,/break-inside:avoid-page/);
});
