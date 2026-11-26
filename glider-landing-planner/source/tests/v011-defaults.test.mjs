import {previouslyPlannedProject} from './previously-planned-fixture.mjs';
import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {createDemo, INITIAL_MAP_ZOOM, INITIAL_LANDING_POINT, makeStraight, makeTurn, ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {effectiveLegSpeedMps} from '../../assets/js/ui/leg-display.js';
import {heightInputText} from '../../assets/js/ui/height-units.js';
import {buildPlan} from '../../assets/js/planning/plan.js';
import {serializeProject, parseProject, loadProject} from '../../assets/js/io/project.js';

function near(a, b, tolerance=1e-8) {
    assert.ok(Math.abs(a-b)<tolerance, `${a} differs from ${b} (tolerance ${tolerance})`);
}

// The main browser module is inert with no application roots. Inspect App's
// real new-plan action without adding production-only test hooks.
globalThis.document={querySelectorAll:()=>[]};
const {App}=await import('../../assets/js/main.js');

test('v0.1.11 new location has the required WGS84 reference coordinate and elevation',()=>{
    const p=createDemo();
    assert.deepEqual(p.landing.point, INITIAL_LANDING_POINT);
    near(p.landing.point.latDeg, 36.17943, 1e-10);
    near(p.landing.point.lonDeg, 137.31361, 1e-10);
    assert.equal(INITIAL_MAP_ZOOM,14);
    near(p.landing.elevationMSLM,714.1464,1e-12);
    assert.equal(heightInputText(p.landing.elevationMSLM,'ft'),'2343');
});

test('v0.1.11 departure altitude is 1200 ft ABOVE reference elevation, not 1200 ft MSL',()=>{
    const p=createDemo();
    near(p.departure.altitudeMSLM - p.landing.elevationMSLM,1200*0.3048,1e-10);
    near(p.departure.altitudeMSLM,1079.9064,1e-10);
    assert.equal(heightInputText(p.departure.altitudeMSLM,'ft'),'3543');
    assert.equal(p.departure.headingDeg,95);
    assert.equal(p.landing.runwayTrackDeg,95);
});

test('new plan preserves both explicit departure coordinates',()=>{
    const p=createDemo();
    near(p.departure.latDeg,36.17914,1e-12);
    near(p.departure.lonDeg,137.31855,1e-12);
});

test('single-model reference performance uses 25 m/s and L/D 34',()=>{
    const p=createDemo();
    assert.equal(p.defaultSpeedMps,25);
    assert.deepEqual(p.performance,{bestLD:34,bestLDAtMps:25});
    near(p.defaultSpeedMps/p.performance.bestLD,25/34,1e-12);
    assert.equal(p.maxBankDeg,60);
    assert.equal(p.wind.fromDeg,270);
    assert.equal(p.wind.speedMps,0);
    assert.equal(p.landing.final.pathAngleDeg,6);
    assert.equal(p.landing.final.groundLengthM,600);
    assert.equal(p.display.magneticDeclinationDeg,0);
});

test('v0.1.11 newly added straight and turn inherit 90 km/h while turn retains bank 30',()=>{
    const p=createDemo();
    for (const leg of [makeStraight(250),makeStraight(),makeTurn()]){
        assert.equal(leg.speedOverrideMps,null);
        near(effectiveLegSpeedMps(p,leg),25);
    }
    assert.equal(makeTurn().bankDeg,30);
});

test('v0.1.11 inherited speeds follow baseline but explicit overrides do not',()=>{
    const p=createDemo();
    const inherited=makeStraight();const explicit=makeTurn();
    explicit.speedOverrideMps=80/3.6;
    p.departureLegs=[inherited,explicit];
    p.defaultSpeedMps=110/3.6;
    near(effectiveLegSpeedMps(p,inherited),110/3.6);
    near(effectiveLegSpeedMps(p,explicit),80/3.6);
    const saved=parseProject(serializeProject(p));
    assert.equal(saved.departureLegs[0].speedOverrideMps,null);
    near(saved.departureLegs[1].speedOverrideMps,80/3.6);
});

test('v0.1.11 changes to L/D and departure point change the plan, without any special-case aircraft',()=>{
    const p=createDemo();
    p.departureLegs=[makeStraight(250)];
    const initial=buildPlan(p);
    p.performance.bestLD=20;
    const changedLd=buildPlan(p);
    assert.notEqual(changedLd.departure.end.h,initial.departure.end.h);
    p.departure.latDeg+=0.002;
    const moved=buildPlan(p);
    assert.notDeepEqual(moved.departure.end,changedLd.departure.end);
});

test('v0.1.11 JSON round trip keeps exact current new-plan values',()=>{
    const p=createDemo();
    const restored=parseProject(serializeProject(p));
    assert.deepEqual(restored,p);
    assert.equal(restored.schemaVersion,'1.3');
    assert.equal(restored.engineVersion,ENGINE_VERSION);
});


test('auto-save ingress preserves a current snapshot',async()=>{
    const project=createDemo();
    project.departure.altitudeMSLM=212.3;
    project.landing.elevationMSLM=33.3;
    project.defaultSpeedMps=24.25;
    project.performance.bestLD=29;
    project.departure.latDeg=34.44;
    const before=structuredClone(project);
    const oldWindow=globalThis.window;
    const oldIndexedDB=globalThis.indexedDB;
    try {
        const db={close(){},transaction(){return {objectStore(){return {openCursor(){
            const req={};
            queueMicrotask(()=>{req.result={value:{project}};req.onsuccess?.();});
            return req;
        }}}};}};
        globalThis.window={indexedDB:{}};
        globalThis.indexedDB={open(){const req={result:db};queueMicrotask(()=>req.onsuccess?.());return req;}};
        const {latestSnapshot}=await import('../../assets/js/io/autosave.js');
        const restored=loadProject(await latestSnapshot());
        assert.deepEqual(restored,before);
        assert.deepEqual(project,before);
    } finally {globalThis.window=oldWindow;globalThis.indexedDB=oldIndexedDB;}
});

test('v0.1.11 reset via New uses exact initial handles, new map view and retains undo',()=>{
    const app=Object.create(App.prototype);
    app.project=previouslyPlannedProject();
    app.project.departure.altitudeMSLM=180;
    app.project.connector.handleStartM=27;
    app.plan=buildPlan(app.project);
    app.history=[];app.future=[];app.selected='old-leg';app.pick='departure';
    app.mobilePanel={kind:'leg',id:'old-leg'};
    app.editRevision=0;
    const mapEvents=[];
    app.map={setPick:(pick)=>mapEvents.push(['pick',pick]),showNewProject:(p)=>mapEvents.push(['new',p.landing.point.latDeg,INITIAL_MAP_ZOOM])};
    app.render=()=>{};app.scheduleSave=()=>{};
    const previous=structuredClone(app.project);
    const originalWindow=globalThis.window;
    globalThis.window={confirm:()=>true};
    try {app.handleAction('new',{});} finally {globalThis.window=originalWindow;}
    assert.equal(app.history.length,1);
    assert.deepEqual(app.history[0],previous);
    near(app.project.departure.altitudeMSLM,1079.9064,1e-10);
    assert.deepEqual(app.project.departureLegs,[]);
    assert.equal(app.project.connector,null);
    assert.equal(app.selected,null);
    assert.equal(app.pick,null);
    assert.equal(app.mobilePanel,null);
    assert.deepEqual(mapEvents,[['pick',null],['new',INITIAL_LANDING_POINT.latDeg,14]]);
    assert.equal(app.plan.connector,null);
    assert.equal(app.editRevision,1);
    assert.equal(app.future.length,0);
});

test('v0.1.11 cancelling New leaves project untouched',()=>{
    const app=Object.create(App.prototype);app.project=createDemo();app.editRevision=9;
    const before=structuredClone(app.project);
    const originalWindow=globalThis.window;globalThis.window={confirm:()=>false};
    try {app.handleAction('new',{});} finally {globalThis.window=originalWindow;}
    assert.deepEqual(app.project,before);
    assert.equal(app.editRevision,9);
});

test('v0.1.11 map constructor centers on the reference with zoom14, not fitBounds',async()=>{
    const previousL=globalThis.L;
    const views=[];let fits=0;
    globalThis.L={
        map(){return {setView(center,zoom){views.push({center,zoom});return this;},fitBounds(){fits++;},on(){return this;},invalidateSize(){}};},
        svg(){return {};},tileLayer(){return {addTo(){return this;}};},layerGroup(){return {addTo(){return this;}};}
    };
    try {
        const {PlannerMap}=await import('../../assets/js/ui/map.js');
        const map=new PlannerMap({innerHTML:''},{onPick(){},onLegLength(){},onFinalLength(){},onHandle(){}});
        assert.deepEqual(views,[{center:[INITIAL_LANDING_POINT.latDeg,INITIAL_LANDING_POINT.lonDeg],zoom:14}]);
        map.showNewProject(createDemo());
        assert.equal(views.length,2);
        assert.equal(fits,0);
    } finally {globalThis.L=previousL;}
});

test('v0.1.11 release and all relative compiled ES modules use the same cache version',()=>{
    const read=(path)=>readFileSync(new URL(path,import.meta.url),'utf8');
    assert.equal(ENGINE_VERSION,'1.0.0');
    assert.match(read('../../glider-landing-planner.php'),/Version:\s*1\.0\.0/);
    assert.match(read('../src/main.ts'),/v1\.0\.0/);
    assert.match(read('../src/ui/print.ts'),/ver=1\.0\.0/);
    assert.match(read('../../README.md'),/v1\.0\.0/);
    assert.match(read('../../README.md'),/https:\/\/mtkbirdman\.com\/glider-landing-planner/);
    assert.match(read('../../preview-local.html'),/assets\/js\/main\.js\?v=1\.0\.0/);
    assert.match(read('../../preview-local.html'),/assets\/css\/app\.css\?v=1\.0\.0/);
    assert.match(read('../../LICENSE'),/GNU GENERAL PUBLIC LICENSE/);
    assert.equal(JSON.parse(read('../package.json')).version,'1.0.0');
    const root=fileURLToPath(new URL('../../assets/js/',import.meta.url));
    function collect(dir){return readdirSync(dir,{withFileTypes:true}).flatMap(entry=>entry.isDirectory()?collect(join(dir,entry.name)):[join(dir,entry.name)]);}
    let modules=0;
    for(const f of collect(root).filter(f=>f.endsWith('.js'))){
        const content=readFileSync(f,'utf8');
        const imports=[...content.matchAll(/\bfrom\s*['"](\.{1,2}\/[^'"]+\.js)(\?v=[^'"]+)?['"]/g)];
        for(const [,mod,version] of imports){
            assert.equal(version,'?v=1.0.0',`${f}: ${mod}`);
            modules++;
        }
    }
    assert.ok(modules>=20,`Expected at least 20 relative module imports, got ${modules}`);
});

test('v0.1.11 user-visible labels do not name an aircraft or initial airfield',()=>{
    const paths=['../src/main.ts','../src/ui/print.ts','../src/ui/map.ts'];
    for(const path of paths){
        const content=readFileSync(new URL(path,import.meta.url),'utf8');
        assert.doesNotMatch(content,/ASK.?21|妻沼滑空場|初期機体|初期滑空場/);
    }
});
