export type GuidanceMode = 'heading'|'ground-track';
export type ModelKind = 'fixed-ld'|'quadratic-drag';
export type Side = 'departure'|'arrival';
export type Direction = 'left'|'right';
export interface GeoPoint {latDeg:number;lonDeg:number}
export interface Point2 {e:number;n:number}
export interface FlightAnchor extends GeoPoint {altitudeMSLM:number;headingDeg:number}
export interface WindConfig {fromDeg:number; speedMps:number}
export interface PerformanceConfig {model:ModelKind;fixedLD:number;bestLD:number;bestLDAtMps:number}
export interface LegBase {id:string;guidance:GuidanceMode;speedOverrideMps:number|null;checkpointLabel?:string;directionDeg:number|null}
export interface StraightLeg extends LegBase {type:'straight';lengthM:number}
export interface TurnLeg extends LegBase {type:'turn';direction:Direction;angleDeg:number;bankDeg:number}
export type FlightLeg = StraightLeg|TurnLeg;
export interface LandingConfig {point:GeoPoint;elevationMSLM:number;runwayTrackDeg:number;final:{pathAngleDeg:number;groundLengthM:number}}
export interface BezierConnector {type:'bezier_connector';handleStartM:number;handleEndM:number;curvaturePolicy:'match-neighbors'}
export interface DisplayConfig {altitudeUnit:'m'|'ft';showDistanceRings:boolean;showBearingLines:boolean;showWindArrow:boolean;showHandles:boolean}
export interface ProjectV1 {schemaVersion:'1.0';engineVersion:string;projectId:string;name:string;defaultSpeedMps:number;maxBankDeg:number;performance:PerformanceConfig;wind:WindConfig;departure:FlightAnchor;landing:LandingConfig;departureLegs:FlightLeg[];arrivalLegs:FlightLeg[];connector:BezierConnector|null;display:DisplayConfig}
export interface Sample {e:number;n:number;h:number;headingDeg:number;trackDeg:number;bankDeg:number;speedMps:number;groundMps:number;tS:number;sM:number; segmentId:string;group:'departure'|'arrival'|'connector'|'final'}
export interface State extends Omit<Sample,'segmentId'|'group'> {}
export interface Notice {code:string;level:'info'|'warning'|'error';message:string;legId?:string}
export interface SegmentResult {leg:FlightLeg;samples:Sample[];start:State;end:State}
export interface PathResult {samples:Sample[];segments:SegmentResult[];start:State;end:State;notices:Notice[]}
export interface ConnectionResult {points:Point2[];samples:Sample[];end:State;lengthM:number;lossM:number;maxBankDeg:number;feasible:boolean;notices:Notice[];handleStartM:number;handleEndM:number}
export interface Checkpoint {label:string;group:'departure'|'connector'|'arrival'|'final';e:number;n:number;predictedM:number;requiredM:number;marginM:number;}
export interface PlanResult {departure:PathResult;arrival:PathResult;connector:ConnectionResult|null;forwardAll:Sample[];requiredAll:Sample[];finalTarget:Sample[];checkpoints:Checkpoint[];notices:Notice[];arrivalStart:State;finalStart:State;landing:Point2;valid:boolean}
