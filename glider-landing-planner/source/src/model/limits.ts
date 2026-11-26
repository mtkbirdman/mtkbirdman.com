/** Persisted connector handle contract (not the local-frame radius). */
export const HANDLE_MIN_M=1, HANDLE_MAX_M=50000;
export function validHandle(value:unknown):value is number {
    return typeof value==='number'&&Number.isFinite(value)&&value>=HANDLE_MIN_M&&value<=HANDLE_MAX_M;
}
export function clampHandle(value:number):number {
    if(!Number.isFinite(value)) throw new Error('Bezierハンドル長が非有限です。');
    return Math.min(HANDLE_MAX_M,Math.max(HANDLE_MIN_M,value));
}
