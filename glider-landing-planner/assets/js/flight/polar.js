import { rad } from '../geometry/geo.js?v=0.1.16';
export class PolarEvaluationError extends Error {
    code;
    constructor(code, message) {
        super(message);
        this.code = code;
        this.name = 'PolarEvaluationError';
    }
}
/** Fritsch–Butland PCHIP slopes: no interpolation beyond supplied speeds. */
export function pchip(points) {
    const n = points.length;
    const x = points.map(p => p.speedMps), y = points.map(p => p.sinkMps);
    const h = Array.from({ length: n - 1 }, (_, i) => x[i + 1] - x[i]);
    const delta = h.map((step, i) => (y[i + 1] - y[i]) / step);
    const slope = Array(n).fill(0);
    for (let i = 1; i < n - 1; i++) {
        const left = delta[i - 1], right = delta[i];
        if (left * right > 0) {
            const w1 = 2 * h[i] + h[i - 1], w2 = h[i] + 2 * h[i - 1];
            slope[i] = (w1 + w2) / (w1 / left + w2 / right);
        }
    }
    // Endpoint derivatives are one-sided and limited to preserve local shape.
    const endpoint = (h0, h1, d0, d1) => {
        let m = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1);
        if (Math.sign(m) !== Math.sign(d0))
            m = 0;
        else if (Math.sign(d0) !== Math.sign(d1) && Math.abs(m) > 3 * Math.abs(d0))
            m = 3 * d0;
        return m;
    };
    slope[0] = endpoint(h[0], h[1], delta[0], delta[1]);
    slope[n - 1] = endpoint(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3]);
    return (v) => {
        if (!Number.isFinite(v) || v < x[0] - 1e-10 || v > x[n - 1] + 1e-10) {
            throw new PolarEvaluationError('POLAR_SPEED_OUT_OF_RANGE', `ポーラー参照速度 ${(v * 3.6).toFixed(1)} km/h はデータ範囲 ${(x[0] * 3.6).toFixed(1)}〜${(x[n - 1] * 3.6).toFixed(1)} km/h の外です。`);
        }
        if (v <= x[0])
            return y[0];
        if (v >= x[n - 1])
            return y[n - 1];
        let low = 0, high = n - 1;
        while (low + 1 < high) {
            const mid = Math.floor((low + high) / 2);
            if (x[mid] <= v)
                low = mid;
            else
                high = mid;
        }
        const t = (v - x[low]) / h[low];
        const t2 = t * t, t3 = t2 * t;
        const result = (2 * t3 - 3 * t2 + 1) * y[low] +
            (t3 - 2 * t2 + t) * h[low] * slope[low] +
            (-2 * t3 + 3 * t2) * y[low + 1] +
            (t3 - t2) * h[low] * slope[low + 1];
        if (!Number.isFinite(result) || result <= 0)
            throw new PolarEvaluationError('POLAR_RESULT_INVALID', 'ポーラーの補間結果が正の有限値になりません。');
        return result;
    };
}
/** One evaluator per plan; the quadratic expression is unchanged from v0.1.11. */
export function createSinkModel(p) {
    if (p.model === 'quadratic-drag') {
        const Vbest = p.quadratic.bestLDAtMps, E = p.quadratic.bestLD;
        const a = 1 / (2 * E * Vbest ** 2), b = Vbest ** 2 / (2 * E);
        return { rate: (v, bank) => a * v ** 3 + b / (v * Math.cos(rad(bank)) ** 2) };
    }
    if (!p.polarTable || p.polarTable.points.length < 4)
        throw new PolarEvaluationError('POLAR_DATA_MISSING', 'ポーラーデータが不足しています。');
    const polar = pchip(p.polarTable.points);
    return { rate: (v, bank) => {
            const n = 1 / Math.cos(rad(bank));
            const result = n ** 1.5 * polar(v / Math.sqrt(n));
            if (!Number.isFinite(result) || result <= 0)
                throw new PolarEvaluationError('POLAR_RESULT_INVALID', '沈下率が正の有限値になりません。');
            return result;
        } };
}
