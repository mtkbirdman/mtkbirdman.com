// Prior releases created a 250 m departure leg and an automatic Bezier.
// Use this fixture when a legacy regression needs those explicit inputs;
// createDemo() now intentionally returns an unplanned project.
import {LocalFrame, headingVector, mul} from '../../assets/js/geometry/geo.js';
import {createDemo,makeStraight} from '../../assets/js/model/defaults.js';
export function previouslyPlannedProject(){
  const p=createDemo();
  p.landing={point:{latDeg:36.211388889,lonDeg:139.418888889},elevationMSLM:89*.3048,runwayTrackDeg:140,final:{pathAngleDeg:6,groundLengthM:450}};
  p.departure={...new LocalFrame(p.landing.point).toGeo(mul(headingVector(320),850)),altitudeMSLM:1289*.3048,headingDeg:140};
  p.departureLegs=[{...makeStraight(250),guidance:'heading'}];
  p.connector={type:'bezier_connector',handleStartM:70,handleEndM:70,curvaturePolicy:'match-neighbors',adjustmentMode:'auto'};
  return p;
}
