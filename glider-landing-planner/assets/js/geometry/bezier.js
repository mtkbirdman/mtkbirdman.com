import { add, sub, mul, magnitude, cross, leftNormal, dot, inLocalRange } from './geo.js?v=0.1.21';
const C = [1, 5, 10, 10, 5, 1];
export function point(cp, t) { let e = 0, n = 0; for (let i = 0; i < 6; i++) {
    const k = C[i] * (1 - t) ** (5 - i) * t ** i;
    e += k * cp[i].e;
    n += k * cp[i].n;
} return { e, n }; }
export function derivative(cp, t) { const dif = cp.slice(1).map((p, i) => mul(sub(p, cp[i]), 5)); let e = 0, n = 0; for (let i = 0; i < 5; i++) {
    const k = [1, 4, 6, 4, 1][i] * (1 - t) ** (4 - i) * t ** i;
    e += k * dif[i].e;
    n += k * dif[i].n;
} return { e, n }; }
export function second(cp, t) { const dif = cp.slice(2).map((p, i) => mul(add(sub(p, mul(cp[i + 1], 2)), cp[i]), 20)); let e = 0, n = 0; for (let i = 0; i < 4; i++) {
    const k = [1, 3, 3, 1][i] * (1 - t) ** (3 - i) * t ** i;
    e += k * dif[i].e;
    n += k * dif[i].n;
} return { e, n }; }
export function curvature(cp, t) { const v = derivative(cp, t), a = second(cp, t); return cross(v, a) / (magnitude(v) ** 3); }
export function controlPoints(A, B, tA, tB, kA, kB, a, b) {
    const p0 = A, p5 = B;
    const p1 = add(A, mul(tA, a));
    const p2 = add(add(A, mul(tA, 2 * a)), mul(leftNormal(tA), (5 / 4) * kA * a * a));
    const p4 = sub(B, mul(tB, b));
    const p3 = add(sub(B, mul(tB, 2 * b)), mul(leftNormal(tB), (5 / 4) * kB * b * b));
    return [p0, p1, p2, p3, p4, p5];
}
const TOL = 1e-6;
export function linesCross(a, b, c, d) {
    const overlap = (a, b, c, d) => Math.max(Math.min(a, b), Math.min(c, d)) <= Math.min(Math.max(a, b), Math.max(c, d)) + TOL;
    if (!overlap(a.e, b.e, c.e, d.e) || !overlap(a.n, b.n, c.n, d.n))
        return false;
    const orient = (p, q, r) => { const v = sub(q, p), z = cross(v, sub(r, p)); return Math.abs(z) <= TOL * magnitude(v) ? 0 : Math.sign(z); };
    return orient(a, b, c) * orient(a, b, d) <= 0 && orient(c, d, a) * orient(c, d, b) <= 0;
}
export function selfIntersects(points) {
    const segments = [];
    for (let i = 0; i < points.length - 1; i++) {
        const a = points[i], b = points[i + 1], v = sub(b, a);
        if (magnitude(v) <= TOL)
            return true;
        if (i + 2 < points.length) {
            const w = sub(points[i + 2], b);
            if (dot(v, w) < 0 && Math.abs(cross(v, w)) <= TOL * Math.max(magnitude(v), magnitude(w)))
                return true;
        }
        segments.push({ i, a, b, minE: Math.min(a.e, b.e), maxE: Math.max(a.e, b.e), minN: Math.min(a.n, b.n), maxN: Math.max(a.n, b.n) });
    }
    segments.sort((a, b) => a.minE - b.minE || a.i - b.i);
    let active = [];
    for (const s of segments) {
        active = active.filter(t => t.maxE + TOL >= s.minE);
        for (const t of active) {
            if (Math.abs(s.i - t.i) <= 1 || s.minN > t.maxN + TOL || t.minN > s.maxN + TOL)
                continue;
            if (linesCross(s.a, s.b, t.a, t.b))
                return true;
        }
        active.push(s);
    }
    return false;
}
/** Convex-hull certificate with bounded subdivision. */
export function curveInRange(cp, depth = 0) {
    if (cp.every(inLocalRange))
        return true;
    if (!inLocalRange(cp[0]) || !inLocalRange(cp.at(-1)) || depth >= 12)
        return false;
    const left = [cp[0]], right = [cp.at(-1)];
    let row = cp;
    while (row.length > 1) {
        row = row.slice(1).map((v, i) => mul(add(v, row[i]), .5));
        left.push(row[0]);
        right.unshift(row.at(-1));
    }
    return curveInRange(left, depth + 1) && curveInRange(right, depth + 1);
}
export const TANGENT_MIN_M = 1e-5;
export const REGULARITY_MAX_DEPTH = 16, REGULARITY_NODE_BUDGET = 20000;
/** Sufficient convex-hull certificate on the original parameter derivative. */
export function curveRegular(cp, maxDepth = REGULARITY_MAX_DEPTH, budget = REGULARITY_NODE_BUDGET) {
    if (cp.length !== 6 || !cp.every(p => Number.isFinite(p.e) && Number.isFinite(p.n)))
        return 'uncertain';
    let remaining = Math.max(0, Math.min(REGULARITY_NODE_BUDGET, Math.floor(budget)));
    const depthLimit = Math.max(0, Math.min(REGULARITY_MAX_DEPTH, Math.floor(maxDepth)));
    const visit = (q, depth) => {
        if (!(remaining > 0))
            return 'uncertain';
        remaining--;
        if (!q.every(p => Number.isFinite(p.e) && Number.isFinite(p.n)))
            return 'uncertain';
        if (magnitude(q[0]) <= TANGENT_MIN_M || magnitude(q.at(-1)) <= TANGENT_MIN_M)
            return 'below-limit';
        const axis = q.reduce((a, b) => add(a, b), { e: 0, n: 0 }), len = magnitude(axis);
        if (Number.isFinite(len) && len > 0 && q.every(v => dot(v, mul(axis, 1 / len)) > TANGENT_MIN_M))
            return 'regular';
        if (!(depth < depthLimit))
            return 'uncertain';
        const left = [q[0]], right = [q.at(-1)];
        let row = q;
        while (row.length > 1) {
            row = row.slice(1).map((v, i) => mul(add(v, row[i]), .5));
            left.push(row[0]);
            right.unshift(row.at(-1));
        }
        const a = visit(left, depth + 1);
        if (a === 'below-limit')
            return a;
        const b = visit(right, depth + 1);
        return b === 'below-limit' ? b : a === 'regular' && b === 'regular' ? 'regular' : 'uncertain';
    };
    return visit(cp.slice(1).map((p, i) => mul(sub(p, cp[i]), 5)), 0);
}
