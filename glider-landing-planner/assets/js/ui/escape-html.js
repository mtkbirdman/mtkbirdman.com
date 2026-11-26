/** Escape user-controlled text before interpolating it into HTML markup. */
export function escapeHtml(value) {
    return String(value ?? '').replace(/[&<>"']/g, char => ({
        '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
    }[char] ?? char));
}
// HTML parsing normalizes CR and NUL. Encode those IDs (and reserved-prefix
// collisions) losslessly before escaping; ordinary IDs retain their DOM value.
const ID_PREFIX = '~glp-id:';
export function htmlLegId(id) {
    const token = id.startsWith(ID_PREFIX) || /[\r\u0000\uD800-\uDFFF]/.test(id)
        ? ID_PREFIX + Array.from({ length: id.length }, (_, i) => id.charCodeAt(i).toString(16).padStart(4, '0')).join('') : id;
    return escapeHtml(token);
}
export function readLegId(token) {
    if (!token?.startsWith(ID_PREFIX))
        return token;
    const hex = token.slice(ID_PREFIX.length);
    if (!/^(?:[0-9a-f]{4})+$/.test(hex))
        return undefined;
    return hex.match(/.{4}/g).map(v => String.fromCharCode(parseInt(v, 16))).join('');
}
