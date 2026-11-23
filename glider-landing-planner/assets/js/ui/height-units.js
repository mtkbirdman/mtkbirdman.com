export const FEET_PER_METER = 1 / 0.3048;
export function heightForDisplay(meters, unit) {
    return unit === 'ft' ? meters * FEET_PER_METER : meters;
}
export function heightFromDisplay(value, unit) {
    return unit === 'ft' ? value * 0.3048 : value;
}
export function heightInputText(meters, unit) {
    // Only round the displayed value; the project continues to store exact meters.
    return String(Number(heightForDisplay(meters, unit).toFixed(2)));
}
