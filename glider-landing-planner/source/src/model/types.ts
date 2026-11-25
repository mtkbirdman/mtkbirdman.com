export type GuidanceMode = 'heading'|'ground-track';
export type ModelKind = 'quadratic-drag'|'polar-table';
export type Side = 'departure'|'arrival';
export type Direction = 'left'|'right';
export interface GeoPoint {latDeg:number;lonDeg:number}
export interface Point2 {e:number;n:number}
export interface FlightAnchor extends GeoPoint {altitudeMSLM:number;headingDeg:number}
export interface WindConfig {fromDeg:number; speedMps:number}
export interface PolarPoint {speedMps:number;sinkMps:number}
export interface PerformanceConfig {model:ModelKind;quadratic:{bestLD:number;bestLDAtMps:number};polarTable:{points:PolarPoint[]}|null}
export interface LegBase {id:string;guidance:GuidanceMode;speedOverrideMps:number|null;checkpointLabel?:string;checkpointKind?:'downwind-start';directionDeg:number|null}
export interface StraightLeg extends LegBase {type:'straight';lengthM:number}
export interface TurnLeg extends LegBase {type:'turn';direction:Direction;angleDeg:number;bankDeg:number;angleInputMode:'angle'|'target';targetDirectionDeg:number|null}
export type FlightLeg = StraightLeg|TurnLeg;
export interface LandingConfig {point:GeoPoint;elevationMSLM:number;runwayTrackDeg:number;final:{pathAngleDeg:number;groundLengthM:number}}
export interface BezierConnector {type:'bezier_connector';handleStartM:number;handleEndM:number;curvaturePolicy:'match-neighbors';adjustmentMode?:'auto'|'manual'}
export interface DisplayConfig {altitudeUnit:'m'|'ft';showDistanceRings:boolean;showBearingLines:boolean;showWindArrow:boolean;showHandles:boolean;showDownwindMargin:boolean;showCheckpointAltitudes:boolean;referenceLineOpacity:number;distanceRingIntervalM:number;bearingLineIntervalDeg:number;magneticDeclinationDeg:number;distanceRingOpacity?:number;bearingLineOpacity?:number}
export interface ProjectV1 {schemaVersion:'1.2';engineVersion:string;projectId:string;name:string;defaultSpeedMps:number;maxBankDeg:number;performance:PerformanceConfig;wind:WindConfig;departure:FlightAnchor;landing:LandingConfig;departureLegs:FlightLeg[];arrivalLegs:FlightLeg[];connector:BezierConnector|null;display:DisplayConfig}
export interface Sample {e:number;n:number;h:number;headingDeg:number;trackDeg:number;bankDeg:number;speedMps:number;groundMps:number;altitudeValid?:boolean;tS:number;sM:number; segmentId:string;group:'departure'|'arrival'|'connector'|'final'}
export interface State extends Omit<Sample,'segmentId'|'group'> {}
export interface Notice {code:string;level:'info'|'warning'|'error';message:string;legId?:string}
export interface SegmentResult {leg:FlightLeg;samples:Sample[];start:State;end:State}
export interface PathResult {samples:Sample[];segments:SegmentResult[];start:State;end:State;notices:Notice[]}
export interface ConnectionResult {points:Point2[];samples:Sample[];end:State;lengthM:number;lossM:number;maxBankDeg:number;feasible:boolean;altitudeValid?:boolean;notices:Notice[];handleStartM:number;handleEndM:number}
export interface Checkpoint {id:string;label:string;kind:'departure-start'|'departure-end'|'arrival-start'|'final-start'|'landing';group:'departure'|'arrival'|'final';e:number;n:number;altitudeM:number;altitudeType:'predicted'|'required'|'target';legId?:string;checkpointKind?:'downwind-start';predictedM?:number;requiredM?:number;marginM?:number;altitudeValid?:boolean;}
export interface PlanResult {departure:PathResult;arrival:PathResult;connector:ConnectionResult|null;forwardAll:Sample[];requiredAll:Sample[];finalTarget:Sample[];checkpoints:Checkpoint[];notices:Notice[];arrivalStart:State;finalStart:State;landing:Point2;valid:boolean;altitudeValid?:boolean}
