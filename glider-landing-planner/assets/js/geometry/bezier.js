import { add, sub, mul, magnitude, cross, leftNormal } from './geo.js?v=0.1.13';
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
export function linesCross(a, b, c, d) {
    function orient(p, q, r) { return cross(sub(q, p), sub(r, p)); }
    const ab1 = orient(a, b, c), ab2 = orient(a, b, d), cd1 = orient(c, d, a), cd2 = orient(c, d, b);
    return ab1 * ab2 < -1e-9 && cd1 * cd2 < -1e-9;
}
export function selfIntersects(points) {
    for (let i = 0; i < points.length - 1; i++)
        for (let j = i + 2; j < points.length - 1; j++) {
            if (i === 0 && j === points.length - 2)
                continue;
            if (linesCross(points[i], points[i + 1], points[j], points[j + 1]))
                return true;
        }
    return false;
}
