import { clampHandle } from '../model/limits.js?v=0.1.20';
import { LocalFrame, headingVector, dot, sub, add, mul, magneticToTrue, inLocalRange } from '../geometry/geo.js?v=0.1.20';
import { windVector } from '../flight/core.js?v=0.1.20';
import { isCompactLayout } from './responsive.js?v=0.1.20';
import { escapeHtml } from './escape-html.js?v=0.1.20';
import { fmtAlt } from './chart.js?v=0.1.20';
import { INITIAL_LANDING_POINT, INITIAL_MAP_ZOOM } from '../model/defaults.js?v=0.1.20';
const marginText = (m, unit) => `${m > 0 ? '+' : ''}${fmtAlt(m, unit)}`;
export class PlannerMap {
    cb;
    map;
    graphics;
    rings;
    controls;
    vectorRenderer;
    picking = null;
    ready = false;
    constructor(element, cb) {
        this.cb = cb;
        if (typeof L === 'undefined') {
            element.innerHTML = '<div class="glp-map-unavailable">地図ライブラリの読込に失敗しました。インターネット接続・CDN制限をご確認ください。</div>';
            return;
        }
        this.map = L.map(element, { zoomControl: true, preferCanvas: false }).setView([INITIAL_LANDING_POINT.latDeg, INITIAL_LANDING_POINT.lonDeg], INITIAL_MAP_ZOOM);
        // SVG remains necessary in Cocoon, and SVG simplification is disabled for smooth sampled arcs.
        this.vectorRenderer = L.svg({ padding: 0.5 });
        L.tileLayer('https://cyberjapandata.gsi.go.jp/xyz/pale/{z}/{x}/{y}.png', { attribution: '&copy; <a href="https://maps.gsi.go.jp/" target="_blank" rel="noopener">国土地理院</a>', maxZoom: 18, minZoom: 3 }).addTo(this.map);
        this.graphics = L.layerGroup().addTo(this.map);
        this.rings = L.layerGroup().addTo(this.map);
        this.controls = L.layerGroup().addTo(this.map);
        this.map.on('click', (ev) => {
            if (this.picking) {
                this.cb.onPick(this.picking, ev.latlng.lat, ev.latlng.lng);
                this.setPick(null);
            }
        });
        this.ready = true;
        setTimeout(() => this.map.invalidateSize(), 250);
    }
    setPick(mode) {
        this.picking = mode;
        if (this.ready)
            this.map.getContainer().style.cursor = mode ? 'crosshair' : '';
    }
    showNewProject(project) {
        if (this.ready)
            this.map.setView([project.landing.point.latDeg, project.landing.point.lonDeg], INITIAL_MAP_ZOOM);
    }
    fit(project, plan) {
        if (!this.ready)
            return;
        const frame = new LocalFrame(project.landing.point);
        const groups = [plan.departure.samples, plan.connector?.points ?? [], plan.arrival.samples, plan.finalTarget];
        const stride = Math.max(1, Math.ceil(groups.reduce((n, g) => n + g.length, 0) / 120));
        const pts = [...(plan.departureInRange ? [project.departure] : []), project.landing.point,
            ...groups.flatMap(g => g.filter((p, i) => (i % stride === 0 || i === g.length - 1) && inLocalRange(p)).map(p => frame.toGeo(p)))];
        this.map.fitBounds(L.latLngBounds(pts.filter(p => Number.isFinite(p.latDeg) && Number.isFinite(p.lonDeg)).map(p => [p.latDeg, p.lonDeg])), { padding: [45, 45], maxZoom: 16 });
    }
    redraw(project, plan, selectedLegId = null) {
        if (!this.ready)
            return;
        this.graphics.clearLayers();
        this.rings.clearLayers();
        this.controls.clearLayers();
        const frame = new LocalFrame(project.landing.point), loc = (p) => { const g = frame.toGeo(p); return [g.latDeg, g.lonDeg]; };
        const draw = (points, color, weight = 3, dashArray, legId) => {
            if (points.length < 2)
                return;
            const coordinates = points.map(loc);
            if (!coordinates.every(p => p.every(Number.isFinite)))
                return;
            const line = L.polyline(coordinates, { color, weight, opacity: .95, dashArray, lineCap: 'round', lineJoin: 'round', smoothFactor: 0, renderer: this.vectorRenderer, interactive: !!legId, bubblingMouseEvents: false });
            if (legId && this.cb.onSelectLeg) {
                // A wider, nearly transparent hit target makes narrow curved legs tappable on phones.
                const click = (ev) => {
                    if (this.picking) {
                        this.cb.onPick(this.picking, ev.latlng.lat, ev.latlng.lng);
                        this.setPick(null);
                    }
                    else
                        this.cb.onSelectLeg?.(legId);
                };
                line.on('click', click);
                if (isCompactLayout())
                    L.polyline(coordinates, { color, weight: 22, opacity: .002, lineCap: 'round', lineJoin: 'round', smoothFactor: 0, renderer: this.vectorRenderer, bubblingMouseEvents: false }).on('click', click).addTo(this.graphics);
            }
            line.addTo(this.graphics);
        };
        // With no actual departure legs, show a non-interactive 250 m heading guide.
        // This line is not a simulated leg and never contributes to distance or height loss.
        if (!project.departureLegs.length && plan.departureInRange) {
            const origin = frame.toLocal(project.departure);
            const course = project.departure.headingDeg;
            const tip = add(origin, mul(headingVector(course), 250));
            if (inLocalRange(tip)) {
                draw([origin, tip], '#2563eb', 2.5, '7,6');
                draw([add(tip, mul(headingVector(course + 150), 22)), tip,
                    add(tip, mul(headingVector(course - 150), 22))], '#2563eb', 2);
            }
        }
        for (const seg of plan.departure.segments)
            draw(seg.samples, '#2563eb', 3.5, seg.completed ? undefined : '8,5', seg.leg.id);
        if (plan.connector)
            draw(plan.connector.points, '#8b5cf6', 3.5, plan.connector.feasible ? undefined : '8,5');
        for (const seg of plan.arrival.segments)
            draw(seg.samples, '#059669', 3.5, seg.completed ? undefined : '8,5', seg.leg.id);
        const end = plan.finalTarget[plan.finalTarget.length - 1];
        if (end) {
            draw([plan.finalStart, end], '#059669', 4, '9,4');
            if (!project.arrivalLegs.length) {
                // Direction arrow is inset from the landing marker so it stays visible.
                const course = project.landing.runwayTrackDeg;
                const arrowTip = add(plan.finalStart, mul(headingVector(course), project.landing.final.groundLengthM * .75));
                draw([add(arrowTip, mul(headingVector(course + 150), 22)), arrowTip,
                    add(arrowTip, mul(headingVector(course - 150), 22))], '#059669', 2);
            }
        }
        const selected = [...plan.departure.segments, ...plan.arrival.segments].find(seg => seg.leg.id === selectedLegId);
        if (selected && selected.samples.length > 1) {
            // White halo + orange overlay remains visible on every basemap without changing the underlying leg colors.
            const pts = selected.samples.map(loc);
            L.polyline(pts, { color: '#fff', dashArray: selected.completed ? undefined : '8,5', weight: 11, opacity: .95, interactive: false, smoothFactor: 0, lineJoin: 'round', lineCap: 'round', renderer: this.vectorRenderer }).addTo(this.graphics);
            L.polyline(pts, { color: '#f59e0b', dashArray: selected.completed ? undefined : '8,5', weight: 5.5, opacity: 1, interactive: false, smoothFactor: 0, lineJoin: 'round', lineCap: 'round', renderer: this.vectorRenderer }).addTo(this.graphics);
        }
        const depIcon = L.divIcon({ className: 'glp-map-point', html: '<span class="glp-point glp-point-blue">離</span>', iconSize: [25, 25], iconAnchor: [12, 12] });
        const landIcon = L.divIcon({ className: 'glp-map-point', html: '<span class="glp-point glp-point-green">着</span>', iconSize: [25, 25], iconAnchor: [12, 12] });
        const compact = isCompactLayout();
        // One choice menu serves the initial departure marker and desktop extension buttons.
        const showAddMenu = (side, point, offerPointEdit = false) => {
            const label = side === 'departure' ? '離脱側' : '着陸側';
            const menu = document.createElement('div');
            menu.className = 'glp-map-add-menu';
            menu.innerHTML = `<strong>${label}へレグを追加</strong>`;
            for (const [type, caption] of [['straight', '＋ 直線'], ['turn', '＋ 旋回']]) {
                const button = document.createElement('button');
                button.type = 'button';
                button.textContent = caption;
                button.addEventListener('click', () => {
                    this.map.closePopup();
                    this.cb.onAddAtEnd?.(side, type);
                });
                menu.appendChild(button);
            }
            if (offerPointEdit) {
                const button = document.createElement('button');
                button.type = 'button';
                button.textContent = '離脱地点・条件を編集';
                button.addEventListener('click', () => {
                    this.map.closePopup();
                    this.cb.onSelectPoint?.('departure');
                });
                menu.appendChild(button);
            }
            L.popup({ className: 'glp-map-add-popup', closeButton: true })
                .setLatLng(loc(point)).setContent(menu).openOn(this.map);
        };
        L.marker([project.departure.latDeg, project.departure.lonDeg], { icon: depIcon, title: project.departureLegs.length ? '離脱地点：クリックして編集' : '離脱地点：クリックしてレグを追加' })
            .on('click', () => {
            if (this.picking)
                return;
            if (!plan.departureInRange || project.departureLegs.length || !this.cb.onAddAtEnd) {
                this.cb.onSelectPoint?.('departure');
            }
            else if (compact) {
                this.cb.onAddAtEnd('departure');
            }
            else {
                showAddMenu('departure', frame.toLocal(project.departure), true);
            }
        }).addTo(this.controls);
        L.marker([project.landing.point.latDeg, project.landing.point.lonDeg], { icon: landIcon, title: '着陸地点・場周設定を開く' })
            .on('click', () => { if (!this.picking)
            this.cb.onSelectPoint?.('landing'); }).addTo(this.controls);
        // Extension buttons appear only after actual legs exist; initially the
        // two point markers and their dashed directional guides remain uncluttered.
        if (this.cb.onAddAtEnd) {
            const size = compact ? 44 : 36;
            for (const [side, point, color] of [['departure', plan.departure.end, 'blue'], ['arrival', plan.arrival.start, 'green']]) {
                if (!(side === 'departure' ? project.departureLegs.length && plan.departure.completed : project.arrivalLegs.length && plan.arrival.completed))
                    continue;
                const label = side === 'departure' ? '離脱側' : '着陸側';
                const icon = L.divIcon({
                    className: 'glp-map-add-point',
                    html: `<span class="glp-map-add-badge glp-map-add-${color}" aria-label="${label}にレグを追加">＋</span>`,
                    iconSize: [size, size], iconAnchor: [size / 2, size / 2]
                });
                L.marker(loc(point), { icon, keyboard: true, zIndexOffset: 1000, title: `${label}へ経路を追加` })
                    .on('click', () => {
                    if (this.picking)
                        return;
                    if (compact)
                        this.cb.onAddAtEnd?.(side);
                    else
                        showAddMenu(side, point);
                }).addTo(this.controls);
            }
        }
        const smallIcon = (name) => L.divIcon({ className: 'glp-map-point', html: '<span class="glp-handle">' + escapeHtml(name) + '</span>', iconSize: [19, 19], iconAnchor: [9, 9] });
        if (plan.connector && plan.connectorReady && plan.connector.points.length && project.display.showHandles) {
            const A = plan.departure.end, B = plan.arrival.start, da = headingVector(A.trackDeg), db = headingVector(B.trackDeg);
            const p1 = add(A, mul(da, project.connector.handleStartM)), p4 = sub(B, mul(db, project.connector.handleEndM));
            for (const [name, p, base, unit, end] of [['A', p1, A, da, 'start'], ['B', p4, B, mul(db, -1), 'end']]) {
                if (!inLocalRange(p))
                    continue;
                L.polyline([loc(base), loc(p)], { color: '#9333ea', dashArray: '3,5', weight: 1.5, opacity: .8, smoothFactor: 0, renderer: this.vectorRenderer }).addTo(this.controls);
                L.marker(loc(p), { draggable: true, icon: smallIcon(name), title: 'Bezier制御ハンドル' }).on('dragend', (ev) => { const ll = ev.target.getLatLng(), q = frame.toLocal({ latDeg: ll.lat, lonDeg: ll.lng }); const length = dot(sub(q, base), unit); this.cb.onHandle(end, Number.isFinite(length) ? clampHandle(length) : length); }).addTo(this.controls);
            }
        }
        // Merge checkpoints sharing one map coordinate, so downwind start is not covered by a redundant connector marker.
        const buckets = new Map();
        for (const cp of plan.checkpoints) {
            // Until the arrival side is planned, only the two anchors and guides are shown.
            if (!project.arrivalLegs.length && !project.connector)
                continue;
            const key = `${Math.round(cp.e / 2)}:${Math.round(cp.n / 2)}`;
            const existing = buckets.get(key) ?? [];
            existing.push(cp);
            buckets.set(key, existing);
        }
        let labelIndex = 0;
        for (const cluster of buckets.values()) {
            const cp = cluster.find(c => c.checkpointKind === 'downwind-start') ?? cluster[0];
            const downwind = cluster.find(c => c.checkpointKind === 'downwind-start');
            const showHeights = project.display.showCheckpointAltitudes;
            const showMargin = Boolean(downwind && project.display.showDownwindMargin);
            const lines = [];
            // Permanent map labels stay short and horizontal; details remain in the hover tooltip and table.
            if (showHeights) {
                lines.push(`<span class="glp-cp-height">${escapeHtml(cp.altitudeValid === false ? '評価不可' : fmtAlt(cp.altitudeM, project.display.altitudeUnit))}</span>`);
            }
            if (showMargin && downwind) {
                const value = downwind.marginM === undefined ? '評価不可' : escapeHtml(marginText(downwind.marginM, project.display.altitudeUnit));
                lines.push(`<strong class="glp-cp-margin">高度余裕 ${value}</strong>`);
            }
            const hover = cluster.map(c => `${c.label}: ${c.altitudeValid === false ? '評価不可' : fmtAlt(c.altitudeM, project.display.altitudeUnit) + ' MSL'}`).join(' / ');
            const marker = L.circleMarker(loc(cp), { radius: 4, fillColor: downwind ? '#b45309' : '#d97706', fillOpacity: 1, color: '#fff', weight: 1, interactive: true, bubblingMouseEvents: false, renderer: this.vectorRenderer });
            marker.bindTooltip(lines.length ? `<div class="glp-map-cp-content">${lines.join('')}</div>` : escapeHtml(hover), lines.length ? { permanent: true, direction: labelIndex % 2 ? 'bottom' : 'top', offset: [0, labelIndex % 2 ? 8 : -5], opacity: 0.96, className: 'glp-checkpoint-label', interactive: false } : { permanent: false });
            marker.addTo(this.graphics);
            labelIndex++;
        }
        if (project.display.showDistanceRings) {
            const step = project.display.distanceRingIntervalM;
            for (let radius = step; radius <= 2000 + 1e-6; radius += step)
                L.circle([project.landing.point.latDeg, project.landing.point.lonDeg], { radius, color: '#7a8797', weight: 1, opacity: project.display.referenceLineOpacity, fill: false, dashArray: '4,5', interactive: false, renderer: this.vectorRenderer }).addTo(this.rings);
        }
        if (project.display.showBearingLines) {
            // Magnetic declination is east-positive: true bearing = magnetic bearing + declination.
            // The physical flight equations, runway course and wind stay true-north referenced.
            for (let magnetic = 0; magnetic < 360; magnetic += project.display.bearingLineIntervalDeg) {
                const trueCourse = magneticToTrue(magnetic, project.display.magneticDeclinationDeg);
                const vec = headingVector(trueCourse), a = loc(plan.landing), b = loc(add(plan.landing, mul(vec, 2000)));
                L.polyline([a, b], { color: '#94a3b8', weight: .8, opacity: project.display.referenceLineOpacity, dashArray: '3,7', smoothFactor: 0, interactive: false, renderer: this.vectorRenderer }).addTo(this.rings);
            }
        }
        if (project.display.showWindArrow && project.wind.speedMps > 0) {
            const v = windVector(project.wind), unit = mul(v, 1 / Math.hypot(v.e, v.n)), center = add(plan.landing, { e: 150, n: 150 }), tail = add(center, mul(unit, -80)), head = add(center, mul(unit, 80));
            L.polyline([loc(tail), loc(head)], { color: '#c2410c', weight: 3, renderer: this.vectorRenderer }).bindTooltip('風: ' + project.wind.fromDeg + '°から ' + project.wind.speedMps + ' m/s').addTo(this.rings);
            // The arrow tip points downwind; use a wide head so it is visible at the initial zoom.
            const back = mul(unit, -45), wing = { e: -unit.n * 32, n: unit.e * 32 };
            L.polyline([loc(add(head, add(back, wing))), loc(head), loc(add(head, sub(back, wing)))], {
                color: '#c2410c', weight: 3, lineCap: 'round', lineJoin: 'round', smoothFactor: 0,
                interactive: false, renderer: this.vectorRenderer
            }).addTo(this.rings);
        }
    }
}
