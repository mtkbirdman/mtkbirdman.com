const between = (value, min, max) => typeof value === 'number' && Number.isFinite(value) && value >= min && value <= max;
function obj(value) {
    if (!value || typeof value !== 'object' || Array.isArray(value)) {
        throw Error('オブジェクトが必要です。');
    }
    return value;
}
function check(valid, message) {
    if (!valid)
        throw Error(message);
}
function validateLegs(legs, ids, label) {
    check(Array.isArray(legs) && legs.length <= 100, label + ' は100レグ以下で指定してください。');
    for (const input of legs) {
        const leg = obj(input);
        const id = leg.id;
        check(typeof id === 'string' && id.length > 0 && id.length <= 80 && !ids.has(id), label + ': レグIDが不正または重複しています。');
        ids.add(id);
        check(['heading', 'ground-track'].includes(String(leg.guidance)), label + ': 飛行方式が不正です。');
        check(leg.speedOverrideMps === null || between(leg.speedOverrideMps, 5, 100), label + ': 速度上書きが不正です。');
        check(leg.directionDeg === null || leg.directionDeg === undefined || between(leg.directionDeg, 0, 360), label + ': 方位指定が不正です。');
        if (leg.checkpointLabel !== undefined) {
            check(typeof leg.checkpointLabel === 'string' && leg.checkpointLabel.length <= 70, label + ': チェックポイント名が長すぎます。');
        }
        if (leg.checkpointKind !== undefined) {
            check(leg.checkpointKind === 'downwind-start', label + ': チェックポイント種別が不正です。');
        }
        if (leg.type === 'straight') {
            check(between(leg.lengthM, 1, 50000), label + ': 直線距離が不正です。');
        }
        else if (leg.type === 'turn') {
            check(['left', 'right'].includes(String(leg.direction)), label + ': 旋回方向が不正です。');
            check(typeof leg.angleDeg === 'number' && Number.isFinite(leg.angleDeg) && leg.angleDeg > 0 && between(leg.bankDeg, 1, 70), label + ': 旋回角・バンク角が不正です。');
            check(leg.angleInputMode === 'angle' || leg.angleInputMode === 'target', label + ': 旋回の入力方式が不正です。');
            check(leg.targetDirectionDeg === null || between(leg.targetDirectionDeg, 0, 360), label + ': 旋回の目標方位が不正です。');
            if (leg.angleInputMode === 'target')
                check(between(leg.targetDirectionDeg, 0, 360), label + ': 目標方位を指定してください。');
        }
        else {
            throw Error(label + ': 未対応のレグ種別です。');
        }
    }
}
/** Validate a complete in-memory project, without modifying it. */
export function validateProject(input) {
    const p = obj(input);
    check(p.schemaVersion === '1.3', '未対応のschemaVersionです（対応: 1.3のみ。旧形式は読込不可）。');
    check(typeof p.engineVersion === 'string' && p.engineVersion.length < 70, 'engineVersionが不正です。');
    check(typeof p.projectId === 'string' && p.projectId.length <= 100, 'projectIdが不正です。');
    check(typeof p.name === 'string' && p.name.length <= 120, '計画名は120文字以内です。');
    check(between(p.defaultSpeedMps, 5, 100), '基本速度は5〜100 m/sの範囲にしてください。');
    check(between(p.maxBankDeg, 5, 70), '最大バンク角は5〜70°の範囲にしてください。');
    const performance = obj(p.performance);
    check(Object.keys(performance).length === 2 && between(performance.bestLD, 2, 100) && between(performance.bestLDAtMps, 5, 100), '性能パラメータが不正です。');
    const wind = obj(p.wind);
    check(between(wind.fromDeg, 0, 360) && between(wind.speedMps, 0, 80), '風向・風速が不正です。');
    const departure = obj(p.departure);
    check(between(departure.latDeg, -90, 90) && between(departure.lonDeg, -180, 180), '離脱地点の緯度経度が不正です。');
    check(between(departure.altitudeMSLM, -500, 15000) && between(departure.headingDeg, 0, 360), '離脱高度・機首方位が不正です。');
    const landing = obj(p.landing);
    const landingPoint = obj(landing.point);
    const final = obj(landing.final);
    check(between(landingPoint.latDeg, -90, 90) && between(landingPoint.lonDeg, -180, 180), '着陸地点の座標が不正です。');
    check(between(landing.elevationMSLM, -500, 15000) && between(landing.runwayTrackDeg, 0, 360), '着陸標高・進入方位が不正です。');
    check(between(final.pathAngleDeg, 0.1, 25) && between(final.groundLengthM, 20, 10000), 'ファイナル角度・距離が不正です。');
    const ids = new Set();
    validateLegs(p.departureLegs, ids, '離脱側レグ');
    validateLegs(p.arrivalLegs, ids, '着陸側レグ');
    if (p.connector !== null) {
        const connector = obj(p.connector);
        check(connector.type === 'bezier_connector' && connector.curvaturePolicy === 'match-neighbors', '接続曲線種別が不正です。');
        check(between(connector.handleStartM, 1, 50000) && between(connector.handleEndM, 1, 50000), 'Bezierハンドル長が不正です。');
        check(['auto', 'manual'].includes(String(connector.adjustmentMode)), 'Bezier調整モードが不正です。');
    }
    const display = obj(p.display);
    check(['m', 'ft'].includes(String(display.altitudeUnit)), '高度単位が不正です。');
    for (const key of ['showDistanceRings', 'showBearingLines', 'showWindArrow', 'showHandles', 'showDownwindMargin', 'showCheckpointAltitudes']) {
        check(typeof display[key] === 'boolean', '表示設定 ' + key + ' が不正です。');
    }
    check(between(display.referenceLineOpacity, 0, 1), '不透明度は0〜1で指定してください。');
    check(between(display.distanceRingIntervalM, 25, 2000), '距離リング間隔は25〜2000mの範囲にしてください。');
    check(typeof display.bearingLineIntervalDeg === 'number' && [15, 30, 45, 60, 90, 120, 180].includes(display.bearingLineIntervalDeg), '方位線間隔は15°刻みかつ360°を均等に分割する値にしてください。');
    check(between(display.magneticDeclinationDeg, -30, 30), '磁気偏角は-30〜+30°の範囲にしてください（東偏が正）。');
    return input;
}
/** Accept only the current schema; never migrate or mutate input. */
export function loadProject(input) {
    return validateProject(structuredClone(obj(input)));
}
export function parseProject(text) {
    if (text.length > 2_000_000)
        throw Error('ファイルが大きすぎます（上限2MB）。');
    return loadProject(JSON.parse(text));
}
export function serializeProject(project) {
    return JSON.stringify(validateProject(project), null, 2);
}
export function downloadProject(project) {
    const blob = new Blob([serializeProject(project)], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = url;
    anchor.download = (project.name.replace(/[^\p{L}\p{N}_-]+/gu, '_') || 'glider_plan') + '.glp.json';
    anchor.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
}
