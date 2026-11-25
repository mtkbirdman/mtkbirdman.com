// Prior releases created a 250 m departure leg and an automatic Bezier.
// Use this fixture when a legacy regression needs those explicit inputs;
// createDemo() now intentionally returns an unplanned project.
import {createDemo,makeStraight} from '../../assets/js/model/defaults.js';
export function previouslyPlannedProject(){
  const p=createDemo();
  p.departureLegs=[makeStraight(250)];
  p.connector={type:'bezier_connector',handleStartM:70,handleEndM:70,curvaturePolicy:'match-neighbors',adjustmentMode:'auto'};
  return p;
}
