import { LocalFrame, headingVector, rad, geographicDistance, LOCAL_RADIUS_M, signed } from '../geometry/geo.js?v=0.1.20';
import { airGround, headingForTrack, simulatePath } from '../flight/core.js?v=0.1.20';
import { evaluateConnector, connectorBoundaryValid } from './connector.js?v=0.1.20';
import { createSinkModel, checkedSink } from '../flight/sink.js?v=0.1.20';
export function buildPlan(project) {
    const frame = new LocalFrame(project.landing.point), lp = frame.toLocal(project.landing.point);
    const departureInRange = geographicDistance(project.departure, project.landing.point) <= LOCAL_RADIUS_M;
    const dep = departureInRange ? frame.toLocal(project.departure) : lp;
    const sinkModel = createSinkModel(project.performance);
    const ctx = { speedMps: project.defaultSpeedMps, wind: project.wind, performance: project.performance, sinkModel };
    const d0 = airGround(project.departure.headingDeg, project.defaultSpeedMps, project.wind);
    const departureStart = { ...dep, h: project.departure.altitudeMSLM, headingDeg: project.departure.headingDeg, trackDeg: d0.track, bankDeg: 0, speedMps: project.defaultSpeedMps, groundMps: d0.speed, tS: 0, sM: 0 };
    const departure = departureInRange ? simulatePath(departureStart, project.departureLegs, ctx, 'departure') : { samples: [], segments: [], start: departureStart, end: departureStart, completed: false, notices: [] };
    const final = project.landing.final, heading = headingVector(project.landing.runwayTrackDeg);
    const fstart = { e: lp.e - heading.e * final.groundLengthM, n: lp.n - heading.n * final.groundLengthM };
    const speedFinal = project.arrivalLegs.at(-1)?.speedOverrideMps ?? project.defaultSpeedMps;
    const solved = headingForTrack(project.landing.runwayTrackDeg, speedFinal, project.wind);
    const finalStart = { ...fstart, h: project.landing.elevationMSLM + final.groundLengthM * Math.tan(rad(final.pathAngleDeg)), headingDeg: solved?.headingDeg ?? NaN, trackDeg: project.landing.runwayTrackDeg, bankDeg: 0, speedMps: speedFinal, groundMps: solved?.groundMps ?? NaN, tS: 0, sM: 0 };
    const arrival = solved ? simulatePath(finalStart, project.arrivalLegs, ctx, 'arrival', true) : { samples: [], segments: [], start: finalStart, end: finalStart, completed: false, notices: [] };
    const notices = [...departure.notices, ...arrival.notices];
    if (!departureInRange)
        notices.push({ code: 'LOCAL_FRAME_OUT_OF_RANGE', level: 'error', message: '離脱地点が局所計算領域（50 km）を超えます。' });
    if (!solved)
        notices.push({ code: project.wind.speedMps < speedFinal ? 'GROUND_SPEED_TOO_LOW' : 'WIND_TRACK_INFEASIBLE', level: 'error', message: 'ファイナルの風補正または対地速度下限を満たしません。' });
    const finalSink = checkedSink(sinkModel, speedFinal, 0, notices);
    const finalBoundaryValid = arrival.completed && Math.abs(signed(arrival.end.trackDeg - finalStart.trackDeg)) <= 0.3 && Math.abs(signed(arrival.end.headingDeg - finalStart.headingDeg)) <= 0.3;
    if (arrival.completed && !finalBoundaryValid)
        notices.push({ code: 'BOUNDARY_STATE_INVALID', level: 'error', message: 'ファイナルとの方位境界が一致しません。' });
    for (const l of [...project.departureLegs, ...project.arrivalLegs])
        if (l.type === 'turn' && l.bankDeg > project.maxBankDeg + 1e-6)
            notices.push({ code: 'BANK_LIMIT_EXCEEDED', level: 'error', message: `レグ ${l.id} の指定バンク角 ${l.bankDeg}° が上限 ${project.maxBankDeg}° を超えています。`, legId: l.id });
    const connectorReady = departure.completed && finalBoundaryValid && !!solved && connectorBoundaryValid(departure.end, arrival.start, project);
    const connection = project.connector && connectorReady ? evaluateConnector(departure.end, arrival.start, project, project.connector.handleStartM, project.connector.handleEndM, false, sinkModel) : null;
    if (connection)
        notices.push(...connection.notices);
    else if (project.connector)
        notices.push({ code: 'BOUNDARY_STATE_INVALID', level: 'error', message: '接続境界を評価できません。' });
    else if (!project.arrivalLegs.length)
        notices.push({ code: 'PLAN_NOT_STARTED', level: 'info', message: '経路は未設計です。着陸側レグまたは場周プリセットを追加すると接続曲線を生成します。' });
    else
        notices.push({ code: 'CONNECTOR_MISSING', level: 'error', message: '離脱側と着陸側が接続されていません。必要高度・到達高度の比較は行いません。' });
    const forwardAll = departure.samples.map(v => ({ ...v }));
    if (connection)
        for (const s of connection.samples.slice(1))
            forwardAll.push(s);
    const joined = !!connection?.feasible && connectorReady;
    const predictedB = connection?.end ?? departure.end;
    const requiredAll = arrival.samples.map(s => ({ ...s, tS: s.tS - arrival.start.tS + (joined ? predictedB.tS : 0), sM: s.sM - arrival.start.sM + (joined ? predictedB.sM : 0) }));
    const heightValid = joined && predictedB.altitudeValid !== false && arrival.start.altitudeValid !== false;
    if (joined)
        for (const s of requiredAll)
            forwardAll.push({ ...s, h: predictedB.h + (s.h - arrival.start.h), altitudeValid: heightValid && s.altitudeValid !== false });
    const endTime = requiredAll.at(-1)?.tS ?? 0, endDistance = requiredAll.at(-1)?.sM ?? 0;
    const finalTarget = [];
    for (let i = 0; i <= 48; i++) {
        const dist = i / 48 * final.groundLengthM;
        const s = { ...finalStart, e: fstart.e + heading.e * dist, n: fstart.n + heading.n * dist, h: project.landing.elevationMSLM + (final.groundLengthM - dist) * Math.tan(rad(final.pathAngleDeg)), tS: solved && finalSink !== null ? endTime + dist / solved.groundMps : NaN, sM: endDistance + dist, segmentId: 'final', group: 'final' };
        finalTarget.push(s);
        if (joined && solved && finalSink !== null)
            forwardAll.push({ ...s, h: predictedB.h + finalStart.h - arrival.start.h - (finalSink ?? 0) * dist / solved.groundMps, altitudeValid: heightValid && finalSink !== null });
    }
    requiredAll.push(...finalTarget);
    notices.push({ code: 'FINAL_BRAKE_UNKNOWN', level: 'info', message: `ファイナル${final.pathAngleDeg}°は幾何学的な目標線です。ブレーキ・フレア・実際の進入成立性は未評価。` });
    const finiteResults = forwardAll.every(s => Number.isFinite(s.tS) && Number.isFinite(s.sM) && Number.isFinite(s.h) && s.altitudeValid !== false);
    if (joined && heightValid && finalSink !== null && !finiteResults && !notices.some(n => n.level === 'error'))
        notices.push({ code: 'PLAN_RESULT_INVALID', level: 'error', message: '計画の必須数値を評価できません。' });
    const valid = joined && heightValid && finalSink !== null && finiteResults && !notices.some(n => n.level === 'error');
    const status = valid ? 'valid' : !project.connector && !project.arrivalLegs.length && !notices.some(n => n.level === 'error') ? 'unplanned' : 'invalid';
    const metrics = { totalDistanceM: valid ? forwardAll.at(-1).sM : null, totalTimeS: valid ? forwardAll.at(-1).tS : null, connectionMarginM: valid ? predictedB.h - arrival.start.h : null };
    // Checkpoint heights are independent physical values: departure forecasts vs. landing-side required heights.
    // Do not place a separate "connector end" checkpoint on top of "downwind start".
    const checkpoints = departureInRange ? [{
            id: 'departure-anchor', label: '離脱地点', kind: 'departure-start', group: 'departure',
            e: departure.start.e, n: departure.start.n, altitudeM: departure.start.h, altitudeType: 'predicted', altitudeValid: departure.start.altitudeValid !== false
        }] : [];
    let departureStraightCount = 0, departureTurnCount = 0;
    for (let i = 0; i < departure.segments.length; i++) {
        const seg = departure.segments[i];
        const sequence = seg.leg.type === 'turn' ? ++departureTurnCount : ++departureStraightCount;
        const generated = seg.leg.type === 'turn' ? `離脱側 第${sequence}旋回 終点` : `離脱側 直線${sequence} 終点`;
        if (!seg.completed)
            continue;
        checkpoints.push({ id: 'departure:' + seg.leg.id + ':end', label: seg.leg.checkpointLabel ?? generated, kind: 'departure-end', group: 'departure',
            e: seg.end.e, n: seg.end.n, altitudeM: seg.end.h, altitudeType: 'predicted', altitudeValid: seg.end.altitudeValid !== false, legId: seg.leg.id });
    }
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
        if (!arrival.segments[i].completed)
            continue;
        const cp = { id: 'arrival:' + leg.id + ':start', label, kind: 'arrival-start', group: 'arrival',
            e: req.e, n: req.n, altitudeM: req.h, altitudeType: 'required', altitudeValid: req.altitudeValid !== false, legId: leg.id,
            ...(isDownwind ? { checkpointKind: 'downwind-start' } : {}) };
        // Margin is defined ONLY at the actual downwind-start checkpoint. Display it only if the connector has passed feasibility.
        if (isDownwind && valid && req.altitudeValid !== false) {
            cp.predictedM = predictedB.h + req.h - arrival.start.h;
            cp.requiredM = req.h;
            cp.marginM = cp.predictedM - req.h;
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
    if (valid && metrics.connectionMarginM !== null && metrics.connectionMarginM < -.5 && !downwind)
        notices.push({ code: 'ALTITUDE_DEFICIT', level: 'warning', message: '接続点のモデル上の高度不足があります。' });
    const rank = { error: 0, warning: 1, info: 2 };
    const distinct = [...new Map(notices.map(n => [JSON.stringify([n.code, n.level, n.legId, n.message]), n])).values()].sort((a, b) => rank[a.level] - rank[b.level]);
    return { departure, arrival, connector: connection, forwardAll, requiredAll, finalTarget, checkpoints, notices: distinct, arrivalStart: arrival.start, finalStart, landing: lp, valid, status, metrics, departureInRange, connectorReady, altitudeValid: heightValid && finalSink !== null };
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
            if (segment.completed && !plan.notices.some(n => n.legId === segment.leg.id && n.level === 'error') && saved?.type === 'turn' && saved.angleInputMode === 'target' && segment.leg.type === 'turn')
                saved.angleDeg = segment.leg.angleDeg;
        }
    }
}
