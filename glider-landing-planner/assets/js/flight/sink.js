import { rad } from '../geometry/geo.js?v=1.0.0';
/** Reuse the simple drag-model coefficients throughout one plan. */
export function createSinkModel(p) {
    const Vbest = p.bestLDAtMps, E = p.bestLD;
    const a = 1 / (2 * E * Vbest ** 2), b = Vbest ** 2 / (2 * E);
    return { parametersValid: Number.isFinite(Vbest) && Vbest > 0 && Number.isFinite(E) && E > 0, rate: (v, bank) => a * v ** 3 + b / (v * Math.cos(rad(bank)) ** 2) };
}
/** Model applicability, separate from the pure mathematical polar. */
export function checkedSink(model, v, bank, notices, legId) {
    const w = model.rate(v, bank), ratio = w / v;
    const code = model.parametersValid === false || !Number.isFinite(v) || v <= 0 || !Number.isFinite(bank) || Math.abs(bank) >= 90 || !Number.isFinite(w) || w <= 0 ? 'SINK_RESULT_INVALID' : ratio > 0.2 ? 'SINK_MODEL_OUT_OF_RANGE' : ratio > 0.1 ? 'SINK_MODEL_APPROXIMATION' : null;
    if (code && !notices.some(n => n.code === code && n.legId === legId))
        notices.push({ code, level: code === 'SINK_MODEL_APPROXIMATION' ? 'warning' : 'error', message: code === 'SINK_MODEL_APPROXIMATION' ? '沈下率比が0.1を超え、浅い滑空近似の精度が低下します。' : code === 'SINK_MODEL_OUT_OF_RANGE' ? '沈下率比が0.2を超え、モデル適用範囲外です。' : '沈下率が正の有限値になりません。', ...(legId ? { legId } : {}) });
    return code && code !== 'SINK_MODEL_APPROXIMATION' ? null : w;
}
