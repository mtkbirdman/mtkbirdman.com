export const ENGINE_VERSION = '0.1.17-test';
// Educational initial coordinates, not surveyed runway or touchdown points.
// See docs/THEORY.md for the initial-condition references and limitations.
export const INITIAL_LANDING_POINT = { latDeg: 36.17943, lonDeg: 137.31361 };
export const INITIAL_MAP_ZOOM = 14;
const INITIAL_LANDING_ELEVATION_MSL_M = 2343 * 0.3048;
const INITIAL_DEPARTURE_HEIGHT_ABOVE_LANDING_M = 1200 * 0.3048;
const INITIAL_DEPARTURE_POINT = { latDeg: 36.17914, lonDeg: 137.31855 };
const INITIAL_DEPARTURE_HEADING_DEG = 95;
const INITIAL_RUNWAY_TRACK_DEG = 95;
// New-plan performance: generic educational parameters referenced to the normal ASK 21.
// This is not an aircraft polar or a safety/performance guarantee.
const INITIAL_REFERENCE_AIRSPEED_MPS = 90 / 3.6;
const INITIAL_BEST_LD = 34;
let nextId = 0;
export function newId() { nextId += 1; return 'leg-' + Date.now().toString(36) + '-' + nextId.toString(36); }
export function makeStraight(lengthM = 300) { return { id: newId(), type: 'straight', lengthM, guidance: 'heading', speedOverrideMps: null, directionDeg: null }; }
export function makeTurn(direction = 'right', angleDeg = 90) { return { id: newId(), type: 'turn', direction, angleDeg, bankDeg: 30, angleInputMode: 'angle', targetDirectionDeg: null, guidance: 'heading', speedOverrideMps: null, directionDeg: null }; }
export function createDemo() {
    const landingPoint = { ...INITIAL_LANDING_POINT };
    const departurePoint = { ...INITIAL_DEPARTURE_POINT };
    return {
        schemaVersion: '1.3', engineVersion: ENGINE_VERSION, projectId: 'demo-' + Date.now(), name: '新規計画',
        defaultSpeedMps: INITIAL_REFERENCE_AIRSPEED_MPS, maxBankDeg: 60,
        performance: { bestLD: INITIAL_BEST_LD, bestLDAtMps: INITIAL_REFERENCE_AIRSPEED_MPS },
        wind: { fromDeg: 270, speedMps: 0 },
        departure: { ...departurePoint, altitudeMSLM: INITIAL_LANDING_ELEVATION_MSL_M + INITIAL_DEPARTURE_HEIGHT_ABOVE_LANDING_M, headingDeg: INITIAL_DEPARTURE_HEADING_DEG },
        landing: { point: landingPoint, elevationMSLM: INITIAL_LANDING_ELEVATION_MSL_M, runwayTrackDeg: INITIAL_RUNWAY_TRACK_DEG, final: { pathAngleDeg: 6, groundLengthM: 600 } },
        departureLegs: [], arrivalLegs: [],
        connector: null,
        display: { altitudeUnit: 'ft', showDistanceRings: true, showBearingLines: true, showWindArrow: true,
            showHandles: false, showDownwindMargin: true, showCheckpointAltitudes: false,
            referenceLineOpacity: 0.5, distanceRingIntervalM: 250, bearingLineIntervalDeg: 45, magneticDeclinationDeg: 0 }
    };
}
export function trafficPreset(side) {
    return [
        { ...makeStraight(700), guidance: 'ground-track', checkpointLabel: 'ダウンウィンド開始', checkpointKind: 'downwind-start' },
        { ...makeTurn(side, 90), guidance: 'ground-track', checkpointLabel: 'ベースターン開始' },
        { ...makeStraight(500), guidance: 'ground-track', checkpointLabel: 'ベースレグ開始' },
        { ...makeTurn(side, 90), guidance: 'ground-track', checkpointLabel: 'ファイナルターン開始' }
    ];
}
