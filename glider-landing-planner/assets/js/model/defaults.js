export const ENGINE_VERSION = '0.1.6-test';
let nextId = 0;
export function newId() { nextId += 1; return 'leg-' + Date.now().toString(36) + '-' + nextId.toString(36); }
export function makeStraight(lengthM = 300) { return { id: newId(), type: 'straight', lengthM, guidance: 'heading', speedOverrideMps: null, directionDeg: null }; }
export function makeTurn(direction = 'right', angleDeg = 90) { return { id: newId(), type: 'turn', direction, angleDeg, bankDeg: 30, guidance: 'heading', speedOverrideMps: null, directionDeg: null }; }
export function createDemo() {
    return {
        schemaVersion: '1.0', engineVersion: ENGINE_VERSION, projectId: 'demo-' + Date.now(), name: '新規計画',
        defaultSpeedMps: 80 / 3.6, maxBankDeg: 60, performance: { model: 'fixed-ld', fixedLD: 30, bestLD: 30, bestLDAtMps: 80 / 3.6 }, wind: { fromDeg: 270, speedMps: 0 },
        departure: { latDeg: 35.409, lonDeg: 136.700, altitudeMSLM: 150, headingDeg: 180 },
        landing: { point: { latDeg: 35.400, lonDeg: 136.700 }, elevationMSLM: 40, runwayTrackDeg: 180, final: { pathAngleDeg: 6, groundLengthM: 450 } },
        departureLegs: [makeStraight(250)], arrivalLegs: [], connector: { type: 'bezier_connector', handleStartM: 70, handleEndM: 70, curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto' },
        display: { altitudeUnit: 'ft', showDistanceRings: true, showBearingLines: true, showWindArrow: true, showHandles: false, showDownwindMargin: true, showCheckpointAltitudes: false, referenceLineOpacity: 0.5, distanceRingIntervalM: 250, bearingLineIntervalDeg: 45, magneticDeclinationDeg: 0 }
    };
}
export function trafficPreset(side) {
    return [
        { ...makeStraight(700), checkpointLabel: 'ダウンウィンド開始', checkpointKind: 'downwind-start' },
        { ...makeTurn(side, 90), checkpointLabel: 'ベースターン開始' },
        { ...makeStraight(500), checkpointLabel: 'ベースレグ開始' },
        { ...makeTurn(side, 90), checkpointLabel: 'ファイナルターン開始' }
    ];
}
