import { LocalFrame, headingVector, dot, sub, add, mul } from '../geometry/geo.js';
import { windVector } from '../flight/core.js';
export class PlannerMap {
    element;
    cb;
    map;
    graphics;
    rings;
    controls;
    vectorRenderer;
    picking = null;
    ready = false;
    constructor(element, cb) {
        this.element = element;
        this.cb = cb;
        if (typeof L === 'undefined') {
            element.innerHTML = '<div class="glp-map-unavailable">地図ライブラリの読込に失敗しました。インターネット接続・CDN制限をご確認ください。</div>';
            return;
        }
        this.map = L.map(element, { zoomControl: true, preferCanvas: false }).setView([35.4, 136.7], 14);
        // Draw vectors using SVG. Canvas paths were not visible on the first Cocoon/WordPress test page.
        // SVG also keeps the path elements inspectable in DevTools.
        this.vectorRenderer = L.svg({ padding: 0.5 });
        L.tileLayer('https://cyberjapandata.gsi.go.jp/xyz/pale/{z}/{x}/{y}.png', { attribution: '&copy; <a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>', maxZoom: 18, minZoom: 3 }).addTo(this.map);
        this.graphics = L.layerGroup().addTo(this.map);
        this.rings = L.layerGroup().addTo(this.map);
        this.controls = L.layerGroup().addTo(this.map);
        this.map.on('click', (ev) => { if (this.picking) {
            this.cb.onPick(this.picking, ev.latlng.lat, ev.latlng.lng);
            this.setPick(null);
        } });
        this.ready = true;
        setTimeout(() => this.map.invalidateSize(), 250);
    }
    setPick(mode) { this.picking = mode; if (this.ready)
        this.map.getContainer().style.cursor = mode ? 'crosshair' : ''; }
    fit(project, plan) { if (!this.ready)
        return; const pts = [project.departure, project.landing.point, ...plan.forwardAll.filter((_, i) => i % Math.max(1, Math.ceil(plan.forwardAll.length / 100)) === 0).map(p => new LocalFrame(project.landing.point).toGeo(p))]; this.map.fitBounds(L.latLngBounds(pts.map(p => [p.latDeg, p.lonDeg])), { padding: [45, 45], maxZoom: 16 }); }
    redraw(project, plan) {
        if (!this.ready)
            return;
        this.graphics.clearLayers();
        this.rings.clearLayers();
        this.controls.clearLayers();
        const f = new LocalFrame(project.landing.point), loc = (p) => { const g = f.toGeo(p); return [g.latDeg, g.lonDeg]; };
        const addPath = (points, color, weight = 3, dashArray) => { if (points.length < 2)
            return; L.polyline(points.map(loc), { color, weight, opacity: .95, dashArray, lineCap: 'round', renderer: this.vectorRenderer }).addTo(this.graphics); };
        addPath(plan.departure.samples, '#2563eb', 3.5);
        if (plan.connector)
            addPath(plan.connector.points, '#8b5cf6', 3.5, plan.connector.feasible ? undefined : '8,5');
        addPath(plan.arrival.samples, '#059669', 3.5);
        const end = plan.finalTarget[plan.finalTarget.length - 1];
        addPath([plan.finalStart, end], '#059669', 4, '9,4');
        const depIcon = L.divIcon({ className: 'glp-map-point', html: '<span class="glp-point glp-point-blue">離</span>', iconSize: [25, 25], iconAnchor: [12, 12] });
        const landIcon = L.divIcon({ className: 'glp-map-point', html: '<span class="glp-point glp-point-green">着</span>', iconSize: [25, 25], iconAnchor: [12, 12] });
        L.marker([project.departure.latDeg, project.departure.lonDeg], { draggable: true, icon: depIcon, title: '離脱地点：ドラッグで変更' }).on('dragend', (e) => { const x = e.target.getLatLng(); this.cb.onPick('departure', x.lat, x.lng); }).addTo(this.controls);
        L.marker([project.landing.point.latDeg, project.landing.point.lonDeg], { draggable: true, icon: landIcon, title: '接地点：ドラッグで変更' }).on('dragend', (e) => { const x = e.target.getLatLng(); this.cb.onPick('landing', x.lat, x.lng); }).addTo(this.controls);
        const smallIcon = (name) => L.divIcon({ className: 'glp-map-point', html: '<span class="glp-handle">' + name + '</span>', iconSize: [19, 19], iconAnchor: [9, 9] });
        for (const side of ['departure', 'arrival']) {
            const path = side === 'departure' ? plan.departure : plan.arrival;
            for (const segment of path.segments) {
                if (segment.leg.type !== 'straight')
                    continue;
                const point = side === 'departure' ? segment.end : segment.start;
                L.marker(loc(point), { draggable: true, icon: smallIcon('↔'), title: 'ドラッグで直線レグ長を変更（進行方向へ投影）' }).on('dragend', (ev) => {
                    const ll = ev.target.getLatLng(), grab = f.toLocal({ latDeg: ll.lat, lonDeg: ll.lng });
                    const origin = side === 'departure' ? segment.start : segment.end;
                    const unit = headingVector(segment.start.trackDeg);
                    const vector = side === 'departure' ? sub(grab, origin) : sub(origin, grab);
                    const newLength = Math.max(1, Math.abs(dot(vector, unit)));
                    this.cb.onLegLength(side, segment.leg.id, newLength);
                }).addTo(this.controls);
            }
        }
        L.marker(loc(plan.finalStart), { draggable: true, icon: smallIcon('F'), title: 'ドラッグでファイナル長を変更' }).on('dragend', (ev) => { const ll = ev.target.getLatLng(), p = f.toLocal({ latDeg: ll.lat, lonDeg: ll.lng }); const direct = headingVector(project.landing.runwayTrackDeg); const v = sub(plan.landing, p); this.cb.onFinalLength(Math.max(20, dot(v, direct))); }).addTo(this.controls);
        if (plan.connector && project.display.showHandles) {
            const A = plan.departure.end, B = plan.arrival.start, da = headingVector(A.trackDeg), db = headingVector(B.trackDeg);
            const p1 = add(A, mul(da, project.connector.handleStartM)), p4 = sub(B, mul(db, project.connector.handleEndM));
            for (const [name, p, base, unit, end] of [['A', p1, A, da, 'start'], ['B', p4, B, mul(db, -1), 'end']]) {
                L.polyline([loc(base), loc(p)], { color: '#9333ea', dashArray: '3,5', weight: 1.5, opacity: .8, renderer: this.vectorRenderer }).addTo(this.controls);
                L.marker(loc(p), { draggable: true, icon: smallIcon(name), title: 'Bezier制御ハンドル' }).on('dragend', (ev) => { const ll = ev.target.getLatLng(), q = f.toLocal({ latDeg: ll.lat, lonDeg: ll.lng }); this.cb.onHandle(end, Math.max(1, dot(sub(q, base), unit))); }).addTo(this.controls);
            }
        }
        for (const cp of plan.checkpoints) {
            L.circleMarker(loc(cp), { radius: 4, fillColor: '#d97706', fillOpacity: 1, color: '#fff', weight: 1, renderer: this.vectorRenderer }).bindTooltip(cp.label + ': ' + (cp.marginM / 0.3048).toFixed(0) + ' ft').addTo(this.graphics);
        }
        if (project.display.showDistanceRings) {
            for (const radius of [250, 500, 1000, 2000])
                L.circle([project.landing.point.latDeg, project.landing.point.lonDeg], { radius, color: '#7a8797', weight: 1, opacity: .5, fill: false, dashArray: '4,5', interactive: false, renderer: this.vectorRenderer }).addTo(this.rings);
        }
        if (project.display.showBearingLines) {
            for (let i = 0; i < 360; i += 45) {
                const vec = headingVector(i), a = loc(plan.landing), b = loc(add(plan.landing, mul(vec, 2000)));
                L.polyline([a, b], { color: '#94a3b8', weight: .8, opacity: .5, dashArray: '3,7', interactive: false, renderer: this.vectorRenderer }).addTo(this.rings);
            }
        }
        if (project.display.showWindArrow && project.wind.speedMps > 0) {
            const v = windVector(project.wind), unit = mul(v, 1 / Math.hypot(v.e, v.n)), center = add(plan.landing, { e: 150, n: 150 }), tail = add(center, mul(unit, -80)), head = add(center, mul(unit, 80));
            L.polyline([loc(tail), loc(head)], { color: '#c2410c', weight: 3, renderer: this.vectorRenderer }).bindTooltip('風: ' + project.wind.fromDeg + '°から ' + project.wind.speedMps + ' m/s').addTo(this.rings);
        }
    }
}
