import type {ProjectV1, FlightLeg, GeoPoint} from './types.js';
import {LocalFrame, headingVector, mul} from '../geometry/geo.js';

export const ENGINE_VERSION = '0.1.13-test';

// New-plan location: Menuma reference point (not a surveyed runway threshold or touchdown point).
// See README.md for the published location and elevation reference.
export const INITIAL_LANDING_POINT: Readonly<GeoPoint> = {latDeg: 36.211388889, lonDeg: 139.418888889};
export const INITIAL_MAP_ZOOM = 14;
const INITIAL_LANDING_ELEVATION_MSL_M = 89 * 0.3048;
const INITIAL_DEPARTURE_HEIGHT_ABOVE_LANDING_M = 1200 * 0.3048;
const INITIAL_DEPARTURE_BEARING_DEG = 320;
const INITIAL_DEPARTURE_DISTANCE_M = 850;
const INITIAL_RUNWAY_TRACK_DEG = 140; // Nominal RW14 reference, not a surveyed bearing.

// New-plan performance: generic educational parameters referenced to the normal ASK 21.
// This is not an aircraft polar or a safety/performance guarantee.
const INITIAL_REFERENCE_AIRSPEED_MPS = 90 / 3.6;
const INITIAL_BEST_LD = 34;

let nextId = 0;
export function newId(): string { nextId += 1; return 'leg-' + Date.now().toString(36) + '-' + nextId.toString(36); }
export function makeStraight(lengthM = 300): FlightLeg { return {id: newId(), type: 'straight', lengthM, guidance: 'heading', speedOverrideMps: null, directionDeg: null}; }
export function makeTurn(direction: 'left' | 'right' = 'right', angleDeg = 90): FlightLeg { return {id: newId(), type: 'turn', direction, angleDeg, bankDeg: 30, angleInputMode: 'angle', targetDirectionDeg: null, guidance: 'heading', speedOverrideMps: null, directionDeg: null}; }

export function createDemo(): ProjectV1 {
    const landingPoint = {...INITIAL_LANDING_POINT};
    const frame = new LocalFrame(landingPoint);
    const departurePoint = frame.toGeo(mul(headingVector(INITIAL_DEPARTURE_BEARING_DEG), INITIAL_DEPARTURE_DISTANCE_M));

    return {
        schemaVersion: '1.2', engineVersion: ENGINE_VERSION, projectId: 'demo-' + Date.now(), name: '新規計画',
        defaultSpeedMps: INITIAL_REFERENCE_AIRSPEED_MPS, maxBankDeg: 60,
        performance: {model: 'quadratic-drag', quadratic: {bestLD: INITIAL_BEST_LD, bestLDAtMps: INITIAL_REFERENCE_AIRSPEED_MPS}, polarTable: null},
        wind: {fromDeg: 270, speedMps: 0},
        departure: {...departurePoint, altitudeMSLM: INITIAL_LANDING_ELEVATION_MSL_M + INITIAL_DEPARTURE_HEIGHT_ABOVE_LANDING_M, headingDeg: 140},
        landing: {point: landingPoint, elevationMSLM: INITIAL_LANDING_ELEVATION_MSL_M, runwayTrackDeg: INITIAL_RUNWAY_TRACK_DEG, final: {pathAngleDeg: 6, groundLengthM: 450}},
        departureLegs: [makeStraight(250)], arrivalLegs: [],
        connector: {type: 'bezier_connector', handleStartM: 70, handleEndM: 70, curvaturePolicy: 'match-neighbors', adjustmentMode: 'auto'},
        display: {altitudeUnit: 'ft', showDistanceRings: true, showBearingLines: true, showWindArrow: true,
            showHandles: false, showDownwindMargin: true, showCheckpointAltitudes: false,
            referenceLineOpacity: 0.5, distanceRingIntervalM: 250, bearingLineIntervalDeg: 45, magneticDeclinationDeg: 0}
    };
}

export function trafficPreset(side:'left'|'right'):FlightLeg[]{return [
 {...makeStraight(700),guidance:'ground-track',checkpointLabel:'ダウンウィンド開始',checkpointKind:'downwind-start'},
 {...makeTurn(side,90),guidance:'ground-track',checkpointLabel:'ベースターン開始'},
 {...makeStraight(500),guidance:'ground-track',checkpointLabel:'ベースレグ開始'},
 {...makeTurn(side,90),guidance:'ground-track',checkpointLabel:'ファイナルターン開始'}
];}
