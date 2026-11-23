// Keep the same breakpoint as the mobile editor CSS (including portrait tablets).
export const COMPACT_LAYOUT_QUERY='(max-width:760px), (max-width:1050px) and (orientation:portrait)';
export function isCompactLayout():boolean {return typeof window!=='undefined'&&typeof window.matchMedia==='function'&&window.matchMedia(COMPACT_LAYOUT_QUERY).matches;}
