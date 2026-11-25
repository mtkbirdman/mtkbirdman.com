import { LocalFrame, headingVector, rad } from '../geometry/geo.js?v=0.1.14';
import { airGround, headingForTrack, simulatePath } from '../flight/core.js?v=0.1.14';
import { evaluateConnector } from './connector.js?v=0.1.14';
import { createSinkModel, PolarEvaluationError } from '../flight/polar.js?v=0.1.14';
export function buildPlan(project) {
    const frame = new LocalFrame(project.landing.point), lp = frame.toLocal(project.landing.point);
    const dep = frame.toLocal(project.departure);
    const sinkModel = createSinkModel(project.performance);
    const ctx = { speedMps: project.defaultSpeedMps, wind: project.wind, performance: project.performance, sinkModel };
    const d0 = airGround(project.departure.headingDeg, project.defaultSpeedMps, project.wind);
    const departureStart = { ...dep, h: project.departure.altitudeMSLM, headingDeg: project.departure.headingDeg, trackDeg: d0.track, bankDeg: 0, speedMps: project.defaultSpeedMps, groundMps: d0.speed, tS: 0, sM: 0 };
    const departure = simulatePath(departureStart, project.departureLegs, ctx, 'departure');
    const final = project.landing.final, heading = headingVector(project.landing.runwayTrackDeg);
    const fstart = { e: lp.e - heading.e * final.groundLengthM, n: lp.n - heading.n * final.groundLengthM };
    const lastArrivalLeg = project.arrivalLegs[project.arrivalLegs.length - 1];
    const lastSpeed = lastArrivalLeg?.speedOverrideMps ?? project.defaultSpeedMps;
    const solved = headingForTrack(project.landing.runwayTrackDeg, lastSpeed, project.wind);
    const targetHeading = solved?.headingDeg ?? project.landing.runwayTrackDeg;
    const gd = airGround(targetHeading, lastSpeed, project.wind);
    const finalStart = { ...fstart, h: project.landing.elevationMSLM + final.groundLengthM * Math.tan(rad(final.pathAngleDeg)), headingDeg: targetHeading, trackDeg: project.landing.runwayTrackDeg, bankDeg: 0, speedMps: lastSpeed, groundMps: solved?.groundMps ?? gd.speed, tS: 0, sM: 0 };
    const arrival = simulatePath(finalStart, project.arrivalLegs, ctx, 'arrival', true);
    const notices = [...departure.notices, ...arrival.notices];
    if (!solved)
        notices.push({ code: 'WIND_TRACK_INFEASIBLE', level: 'error', message: 'ファイナルの対地進入方向を維持できる速度・風の組み合わせではありません。' });
    const connection = project.connector ? evaluateConnector(departure.end, arrival.start, project, project.connector.handleStartM, project.connector.handleEndM, false, sinkModel) : null;
    if (connection)
        notices.push(...connection.notices);
    else
        notices.push({ code: 'CONNECTOR_MISSING', level: 'error', message: '離脱側と着陸側が接続されていません。必要高度・到達高度の比較は行いません。' });
    const forwardAll = departure.samples.map(v => ({ ...v }));
    let predictedB = { ...departure.end };
    if (connection) {
        forwardAll.push(...connection.samples.slice(1).map(v => ({ ...v })));
        predictedB = { ...connection.end };
    }
    const predArrivalStart = { ...arrival.start, h: predictedB.h, tS: predictedB.tS, sM: predictedB.sM, ...(predictedB.altitudeValid === false ? { altitudeValid: false } : {}) };
    // Use exactly the turns resolved during reverse planning (including their revolution counts).
    const predictedArrival = simulatePath(predArrivalStart, arrival.segments.map(segment => segment.leg.type === 'turn' ? { ...segment.leg, angleInputMode: 'angle' } : segment.leg), ctx, 'arrival', false);
    notices.push(...predictedArrival.notices.filter(n => n.level === 'error'));
    forwardAll.push(...predictedArrival.samples.slice(1).map(v => ({ ...v })));
    const points = 48;
    const finalTarget = [];
    const finalPredicted = [];
    const speedFinal = project.defaultSpeedMps;
    const finalSol = headingForTrack(project.landing.runwayTrackDeg, speedFinal, project.wind);
    const vg = finalSol?.groundMps ?? airGround(project.landing.runwayTrackDeg, speedFinal, project.wind).speed;
    let finalSink = null;
    try {
        finalSink = sinkModel.rate(speedFinal, 0);
    }
    catch (error) {
        if (!(error instanceof PolarEvaluationError))
            throw error;
        notices.push({ code: error.code, level: 'error', message: error.message });
    }
    const finalLoss = vg > 0.001 && finalSink !== null ? final.groundLengthM * finalSink / vg : 0;
    const predictedFinalValid = predictedArrival.end.altitudeValid !== false && finalSink !== null;
    for (let i = 0; i <= points; i++) {
        const u = i / points, dist = u * final.groundLengthM;
        const e = fstart.e + heading.e * dist, n = fstart.n + heading.n * dist;
        const offsetS = predictedArrival.end.sM + dist, offsetTime = predictedArrival.end.tS + (vg > 0 ? dist / vg : 0);
        finalTarget.push({ e, n, h: project.landing.elevationMSLM + (final.groundLengthM - dist) * Math.tan(rad(final.pathAngleDeg)), headingDeg: finalSol?.headingDeg ?? targetHeading, trackDeg: project.landing.runwayTrackDeg, bankDeg: 0, speedMps: speedFinal, groundMps: vg, tS: offsetTime, sM: offsetS, segmentId: 'final', group: 'final' });
        finalPredicted.push({ e, n, h: predictedArrival.end.h - (dist / final.groundLengthM) * finalLoss, headingDeg: finalSol?.headingDeg ?? targetHeading, trackDeg: project.landing.runwayTrackDeg, bankDeg: 0, speedMps: speedFinal, groundMps: vg, tS: offsetTime, sM: offsetS, segmentId: 'final', group: 'final', ...(predictedFinalValid ? {} : { altitudeValid: false }) });
    }
    forwardAll.push(...finalPredicted.slice(1));
    const requiredAll = [];
    const offsetS = predictedB.sM - arrival.start.sM;
    for (const s of arrival.samples) {
        requiredAll.push({ ...s, sM: s.sM + offsetS });
    }
    const arrivalNeedEndS = arrival.end.sM + offsetS;
    for (const s of finalTarget) {
        requiredAll.push({ ...s, sM: arrivalNeedEndS + (s.sM - predictedArrival.end.sM) });
    }
    // Checkpoint heights are independent physical values: departure forecasts vs. landing-side required heights.
    // Do not place a separate "connector end" checkpoint on top of "downwind start".
    const checkpoints = [{
            id: 'departure-anchor', label: '離脱地点', kind: 'departure-start', group: 'departure',
            e: departure.start.e, n: departure.start.n, altitudeM: departure.start.h, altitudeType: 'predicted', altitudeValid: departure.start.altitudeValid !== false
        }];
    let departureStraightCount = 0, departureTurnCount = 0;
    for (let i = 0; i < departure.segments.length; i++) {
        const seg = departure.segments[i];
        const sequence = seg.leg.type === 'turn' ? ++departureTurnCount : ++departureStraightCount;
        const generated = seg.leg.type === 'turn' ? `離脱側 第${sequence}旋回 終点` : `離脱側 直線${sequence} 終点`;
        checkpoints.push({ id: 'departure:' + seg.leg.id + ':end', label: seg.leg.checkpointLabel ?? generated, kind: 'departure-end', group: 'departure',
            e: seg.end.e, n: seg.end.n, altitudeM: seg.end.h, altitudeType: 'predicted', altitudeValid: seg.end.altitudeValid !== false, legId: seg.leg.id });
    }
    const marginAvailable = Boolean(connection?.feasible) && departure.end.altitudeValid !== false && arrival.start.altitudeValid !== false && !notices.some(n => n.level === 'error');
    let arrivalStraightCount = 0, arrivalTurnCount = 0;
    for (let i = 0; i < arrival.segments.length; i++) {
        const req = arrival.segments[i].start;
        const leg = arrival.segments[i].leg;
        const isDownwind = leg.checkpointKind === 'downwind-start';
        // Preset and user-defined checkpoint names are retained verbatim. Only unnamed legs get numbered labels.
        if (!leg.checkpointLabel) {
            if (leg.type === 'turn')
                arrivalTurnCount++;
            else
                arrivalStraightCount++;
        }
        const generated = leg.type === 'turn' ? `着陸側 第${arrivalTurnCount}旋回 開始` : `着陸側 直線${arrivalStraightCount} 開始`;
        const label = leg.checkpointLabel ?? (isDownwind ? 'ダウンウィンド開始' : generated);
        const cp = { id: 'arrival:' + leg.id + ':start', label, kind: 'arrival-start', group: 'arrival',
            e: req.e, n: req.n, altitudeM: req.h, altitudeType: 'required', altitudeValid: req.altitudeValid !== false, legId: leg.id,
            ...(isDownwind ? { checkpointKind: 'downwind-start' } : {}) };
        // Margin is defined ONLY at the actual downwind-start checkpoint. Display it only if the connector has passed feasibility.
        const pred = predictedArrival.segments[i]?.start;
        if (isDownwind && marginAvailable && pred && pred.altitudeValid !== false && req.altitudeValid !== false) {
            cp.predictedM = pred.h;
            cp.requiredM = req.h;
            cp.marginM = pred.h - req.h;
        }
        checkpoints.push(cp);
    }
    checkpoints.push({ id: 'final-start', label: 'ファイナル開始', kind: 'final-start', group: 'final',
        e: fstart.e, n: fstart.n, altitudeM: finalStart.h, altitudeType: 'target' });
    checkpoints.push({ id: 'landing-target', label: '進入目標点', kind: 'landing', group: 'final',
        e: lp.e, n: lp.n, altitudeM: project.landing.elevationMSLM, altitudeType: 'target' });
    // Two checkpoints can have identical locations (e.g., a zero-length segment). Labels are not duplicated by their name.
    const downwind = checkpoints.find(c => c.checkpointKind === 'downwind-start');
    if (downwind && downwind.marginM !== undefined && downwind.marginM < -0.5)
        notices.push({ code: 'ALTITUDE_DEFICIT', level: 'warning', message: `ダウンウィンド開始: モデル上の高度不足 ${Math.abs(downwind.marginM / 0.3048).toFixed(0)} ft` });
    if (downwind && downwind.marginM === undefined)
        notices.push({ code: 'MARGIN_NOT_EVALUABLE', level: 'info', message: 'ダウンウィンド開始の高度余裕は、有効な接続経路がないため評価できません。' });
    if (!connection) {
        forwardAll.splice(departure.samples.length);
    }
    notices.push({ code: 'FINAL_BRAKE_UNKNOWN', level: 'info', message: `ファイナル${final.pathAngleDeg}°は幾何学的な目標線です。ブレーキ・フレア・実際の進入成立性は未評価。` });
    for (const l of [...project.departureLegs, ...project.arrivalLegs]) {
        if (l.type === 'turn' && l.bankDeg > project.maxBankDeg) {
            notices.push({ code: 'BANK_LIMIT_EXCEEDED', level: 'error', message: `レグ ${l.id} の指定バンク角 ${l.bankDeg}° が上限 ${project.maxBankDeg}° を超えています。`, legId: l.id });
        }
    }
    if (connection && Math.abs(connection.end.headingDeg - arrival.start.headingDeg) > 0.3) {
        notices.push({ code: 'BOUNDARY_DISCONTINUITY', level: 'warning', message: 'Bezierと着陸側経路の機首方位が不連続になる場合があります。' });
    }
    const valid = !notices.some(n => n.level === 'error');
    return { departure, arrival, connector: connection, forwardAll, requiredAll, finalTarget, checkpoints, notices, arrivalStart: arrival.start, finalStart, landing: lp, valid, altitudeValid: departure.end.altitudeValid !== false && arrival.start.altitudeValid !== false && connection?.altitudeValid !== false && predictedFinalValid };
}
/** Save the latest resolved angle as the reference for choosing revolutions on the next edit.
 * Never alter angle-driven turns or the target bearing. */
export function rememberResolvedTurnAngles(project, plan) {
    for (const [legs, segments] of [
        [project.departureLegs, plan.departure.segments],
        [project.arrivalLegs, plan.arrival.segments]
    ]) {
        for (const segment of segments) {
            const saved = legs.find(leg => leg.id === segment.leg.id);
            if (saved?.type === 'turn' && saved.angleInputMode === 'target' && segment.leg.type === 'turn')
                saved.angleDeg = segment.leg.angleDeg;
        }
    }
}
