export type HeightUnit='ft'|'m';
export const FEET_PER_METER=1/0.3048;
export function heightForDisplay(meters:number,unit:HeightUnit):number {
 return unit==='ft'?meters*FEET_PER_METER:meters;
}
export function heightFromDisplay(value:number,unit:HeightUnit):number {
 return unit==='ft'?value*0.3048:value;
}
export function heightInputText(meters:number,unit:HeightUnit):string {
 // Only round the displayed value; the project continues to store exact meters.
 return String(Number(heightForDisplay(meters,unit).toFixed(2)));
}
