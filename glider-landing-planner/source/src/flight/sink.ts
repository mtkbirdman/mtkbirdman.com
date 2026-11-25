import type { PerformanceConfig } from '../model/types.js';
import { rad } from '../geometry/geo.js';

export interface SinkModel {
    rate(speedMps: number, bankDeg: number): number;
}

/** Reuse the simple drag-model coefficients throughout one plan. */
export function createSinkModel(p: PerformanceConfig): SinkModel {
    const Vbest = p.bestLDAtMps, E = p.bestLD;
    const a = 1 / (2 * E * Vbest ** 2), b = Vbest ** 2 / (2 * E);
    return { rate: (v, bank) => a * v ** 3 + b / (v * Math.cos(rad(bank)) ** 2) };
}
