const a = 6378137;
const f = 1 / 298.257223563;
const e2 = f * (2 - f);
const d = Math.PI / 180;
export const rad = (deg) => deg * d;
export const deg = (r) => r / d;
export const norm = (n) => ((n % 360) + 360) % 360;
export const signed = (n) => ((n + 180) % 360 + 360) % 360 - 180;
export const headingVector = (angleDeg) => ({ e: Math.sin(rad(angleDeg)), n: Math.cos(rad(angleDeg)) });
export const courseOf = (e, n) => norm(deg(Math.atan2(e, n)));
export const dot = (a, b) => a.e * b.e + a.n * b.n;
export const cross = (a, b) => a.e * b.n - a.n * b.e;
export const magnitude = (v) => Math.hypot(v.e, v.n);
export const add = (a, b) => ({ e: a.e + b.e, n: a.n + b.n });
export const sub = (a, b) => ({ e: a.e - b.e, n: a.n - b.n });
export const mul = (a, s) => ({ e: a.e * s, n: a.n * s });
export const leftNormal = (v) => ({ e: -v.n, n: v.e });
function xyz(geo, h = 0) { const lat = rad(geo.latDeg), lon = rad(geo.lonDeg); const n = a / Math.sqrt(1 - e2 * Math.sin(lat) ** 2); return { x: (n + h) * Math.cos(lat) * Math.cos(lon), y: (n + h) * Math.cos(lat) * Math.sin(lon), z: (n * (1 - e2) + h) * Math.sin(lat) }; }
export class LocalFrame {
    origin;
    originXYZ;
    slat;
    clat;
    slon;
    clon;
    constructor(origin) {
        this.origin = origin;
        this.originXYZ = xyz(origin);
        this.slat = Math.sin(rad(origin.latDeg));
        this.clat = Math.cos(rad(origin.latDeg));
        this.slon = Math.sin(rad(origin.lonDeg));
        this.clon = Math.cos(rad(origin.lonDeg));
    }
    toLocal(geo) { const p = xyz(geo), dx = p.x - this.originXYZ.x, dy = p.y - this.originXYZ.y, dz = p.z - this.originXYZ.z; return { e: -this.slon * dx + this.clon * dy, n: -this.slat * this.clon * dx - this.slat * this.slon * dy + this.clat * dz }; }
    toGeo(pt) { const { e, n } = pt; const x = this.originXYZ.x - this.slon * e - this.slat * this.clon * n, y = this.originXYZ.y + this.clon * e - this.slat * this.slon * n, z = this.originXYZ.z + this.clat * n; const lon = Math.atan2(y, x), p = Math.hypot(x, y); if (p < 1e-9)
        return { latDeg: z >= 0 ? 90 : -90, lonDeg: this.origin.lonDeg }; let lat = Math.atan2(z, p * (1 - e2)); for (let k = 0; k < 8; k++) {
        const n0 = a / Math.sqrt(1 - e2 * Math.sin(lat) ** 2);
        const h = p / Math.cos(lat) - n0;
        lat = Math.atan2(z, p * (1 - e2 * n0 / (n0 + h)));
    } return { latDeg: deg(lat), lonDeg: deg(lon) }; }
}
/** East-positive magnetic declination; flight-state headings are always true. */
export function magneticToTrue(magneticDeg, eastDeclinationDeg) { return norm(magneticDeg + eastDeclinationDeg); }
export const LOCAL_RADIUS_M = 50000;
export const inLocalRange = (p) => magnitude(p) <= LOCAL_RADIUS_M;
export function geographicDistance(a, b) {
    const unit = (p) => [Math.cos(rad(p.latDeg)) * Math.cos(rad(p.lonDeg)), Math.cos(rad(p.latDeg)) * Math.sin(rad(p.lonDeg)), Math.sin(rad(p.latDeg))];
    const u = unit(a), v = unit(b);
    return 6371008.8 * Math.atan2(Math.hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]), u.reduce((s, x, i) => s + x * v[i], 0));
}
