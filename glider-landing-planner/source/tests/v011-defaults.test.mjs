import test from 'node:test';
import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {join} from 'node:path';
import {fileURLToPath} from 'node:url';
import {createDemo, INITIAL_MAP_ZOOM, INITIAL_LANDING_POINT, makeStraight, makeTurn, ENGINE_VERSION} from '../../assets/js/model/defaults.js';
import {LocalFrame, headingVector, courseOf} from '../../assets/js/geometry/geo.js';
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
    near(p.landing.point.latDeg, 36.211388889, 1e-10);
    near(p.landing.point.lonDeg, 139.418888889, 1e-10);
    assert.equal(INITIAL_MAP_ZOOM,14);
    near(p.landing.elevationMSLM,27.1272,1e-12);
    assert.equal(heightInputText(p.landing.elevationMSLM,'ft'),'89');
});

test('v0.1.11 departure altitude is 1200 ft ABOVE reference elevation, not 1200 ft MSL',()=>{
    const p=createDemo();
    near(p.departure.altitudeMSLM - p.landing.elevationMSLM,1200*0.3048,1e-10);
    near(p.departure.altitudeMSLM,392.8872,1e-10);
    assert.equal(heightInputText(p.departure.altitudeMSLM,'ft'),'1289');
    assert.equal(p.departure.headingDeg,140);
    assert.equal(p.landing.runwayTrackDeg,140);
});

test('v0.1.11 departure location is 850 m / true bearing 320 from landing',()=>{
    const p=createDemo();
    const local=new LocalFrame(p.landing.point).toLocal(p.departure);
    near(Math.hypot(local.e,local.n),850,0.0001);
    near(courseOf(local.e,local.n),320,0.0001);
    const expected=headingVector(320);
    near(local.e/850,expected.e,1e-6);
});

test('v0.1.12 model B reference performance uses 25 m/s and L/D 34',()=>{
    const p=createDemo();
    assert.equal(p.defaultSpeedMps,25);
    assert.deepEqual(p.performance,{model:'quadratic-drag',quadratic:{bestLD:34,bestLDAtMps:25},polarTable:null});
    near(p.defaultSpeedMps/p.performance.quadratic.bestLD,25/34,1e-12);
    assert.equal(p.maxBankDeg,60);
    assert.equal(p.wind.fromDeg,270);
    assert.equal(p.wind.speedMps,0);
    assert.equal(p.landing.final.pathAngleDeg,6);
    assert.equal(p.landing.final.groundLengthM,450);
    assert.equal(p.display.magneticDeclinationDeg,0);
});

test('v0.1.11 newly added straight and turn inherit 90 km/h while turn retains bank 30',()=>{
    const p=createDemo();
    for (const leg of [p.departureLegs[0],makeStraight(),makeTurn()]){
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
    const initial=buildPlan(p);
    p.performance.quadratic.bestLD=20;
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
    assert.equal(restored.schemaVersion,'1.2');
    assert.equal(restored.engineVersion,ENGINE_VERSION);
});

test('v0.1.11 import of a v0.1.10 record preserves all saved fields',()=>{
    // Simulated existing record using the exact v0.1.10 schema and defaults.
    const existing=createDemo();
    existing.engineVersion='0.1.10-test';
    existing.defaultSpeedMps=80/3.6;
    existing.schemaVersion='1.0';
    existing.performance={model:'quadratic-drag',fixedLD:30,bestLD:30,bestLDAtMps:80/3.6};
    existing.departure={latDeg:35.409,lonDeg:136.7,altitudeMSLM:150,headingDeg:180};
    existing.landing={point:{latDeg:35.4,lonDeg:136.7},elevationMSLM:40,runwayTrackDeg:180,final:{pathAngleDeg:6,groundLengthM:450}};
    existing.departureLegs[0].speedOverrideMps=23;
    const before=structuredClone(existing);
    const loaded=parseProject(JSON.stringify(existing));
    assert.equal(loaded.schemaVersion,'1.2');
    assert.deepEqual(loaded.performance,{model:'quadratic-drag',quadratic:{bestLD:30,bestLDAtMps:80/3.6},polarTable:null});
    assert.deepEqual(loaded.departure,before.departure);
    assert.deepEqual(loaded.landing,before.landing);
    assert.deepEqual(existing,before);
    assert.deepEqual(loadProject(existing),loaded);
});

test('v0.1.11 auto-save ingress keeps a v0.1.10 snapshot unchanged',async()=>{
    const project=createDemo();project.engineVersion='0.1.10-test';
    project.departure.altitudeMSLM=212.3;
    project.landing.elevationMSLM=33.3;
    project.defaultSpeedMps=24.25;
    project.performance.quadratic.bestLD=29;
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
    app.project=createDemo();
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
    near(app.project.departure.altitudeMSLM,392.8872,1e-10);
    assert.equal(app.project.connector.handleStartM,70);
    assert.equal(app.project.connector.handleEndM,70);
    assert.equal(app.selected,null);
    assert.equal(app.pick,null);
    assert.equal(app.mobilePanel,null);
    assert.deepEqual(mapEvents,[['pick',null],['new',INITIAL_LANDING_POINT.latDeg,14]]);
    assert.deepEqual(app.plan.connector.handleStartM,buildPlan(createDemo()).connector.handleStartM);
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
    assert.equal(ENGINE_VERSION,'0.1.13-test');
    assert.match(read('../../glider-landing-planner.php'),/Version:\s*0\.1\.13/);
    assert.match(read('../src/main.ts'),/テスト版 v0\.1\.13/);
    assert.match(read('../src/ui/print.ts'),/ver=0\.1\.13/);
    assert.match(read('../../README.md'),/v0\.1\.13/);
    assert.equal(JSON.parse(read('../package.json')).version,'0.1.13');
    const root=fileURLToPath(new URL('../../assets/js/',import.meta.url));
    function collect(dir){return readdirSync(dir,{withFileTypes:true}).flatMap(entry=>entry.isDirectory()?collect(join(dir,entry.name)):[join(dir,entry.name)]);}
    let modules=0;
    for(const f of collect(root).filter(f=>f.endsWith('.js'))){
        const content=readFileSync(f,'utf8');
        const imports=[...content.matchAll(/\bfrom\s*['"](\.{1,2}\/[^'"]+\.js)(\?v=[^'"]+)?['"]/g)];
        for(const [,mod,version] of imports){
            assert.equal(version,'?v=0.1.13',`${f}: ${mod}`);
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
