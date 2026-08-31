# ウインチ曳航シミュレーターの実装仕様

## 1. 文書情報

- 仕様バージョン: **v1.0**
- 対象実装: **winch_launch_simulator v1.0**
- 対象: G103A 縦運動ウインチ発航
- 基準座標: 地面固定 x 前方、z 下方。機体固定 x 前方、z 下方。
- canonical launch distance: **1100 m**

v1.0 はv2.11のClub pilot SOP、3.5 s launch-power ramp、1100 m launch setup、Vw protection、terminal fade、9-state base ODEを維持する。v1.0ではFEMの物理式・reeling/remesh ruleを変更せず、WLO07索空力・内部軸力・索地面接触の配列演算化とRK4中間段のdiagnostics省略で実装を整理し、未完だったfull-launch dt/2・dt/4および20/22/40-node V&Vを完了する。

## 2. モデル階層

同一の地上静止状態から、次の3段階で物理忠実度を上げる。

1. **Prescribed tension**
   - 機体3自由度 + 車輪接地
   - 索・ウインチ動力学は解かず、張力履歴を外部規定
2. **Dynamic winch + Secant cable**
   - ウインチ回転運動、drive torque、制御則
   - 索は等価ばね・ダンパ
   - Standard SOP / Club SOP
3. **Dynamic winch + FEM cable**
   - ウインチ回転運動
   - 分布弾性索、索重量、索空力
   - cable-ground contact / friction
   - reel-in / accepted-step node removal
   - Standard SOP / Club SOP

Standard / Club は fidelity level ではなく、同じ plant に対する操縦・ウインチ操作の差として扱う。

## 3. 解析開始点

canonical case は `launch_setup.mode = "ground_start"` とし、t=0 は次を満たす pre-launch state とする。

- 機体は地上静止
- drum rate = 0
- initial pretension = `2500 N`
- 車輪 ground reaction を含む機体鉛直・pitch平衡
- FEM は重力・地面接触を含む静的平衡形状
- Secant / FEM の material natural length = **1100 m**
- 必要な pre-launch reel-in は initializer が派生量として求める

旧 `h_cg_liftoff_m`、25 m/s airborne start、FEM直線配置、人工的な短縮自然長、非零 initial drum rate は canonical path では使用しない。

## 4. G103A geometry

代表CG基準の body-fixed 座標を使用する。

- belly hook: `(x,z) = (+0.535,+0.366) m`
- nosewheel center: `(+2.029,+0.508) m`
- mainwheel center: `(-0.169,+0.366) m`
- tailwheel center x: `-4.537 m`
- tailwheel center z: `+0.061 m`

タイヤ半径 engineering baseline:

- nose: `0.130 m`
- main: `0.220 m`
- tail: `0.105 m`

hook と nose/main wheel center はユーザー計測値をcanonical geometryとする。tailwheel z、タイヤ半径、loaded radius / load-deflectionは今後の実測更新対象である。

## 5. Aircraft ground contact

flat ground `z = ground.z_m` のみを扱う。wheel center のground-zと半径から signed depth を計算する。

\[
\delta_i=\max(0,z_{c,i}+R_i-z_g)
\]

法線反力:

\[
N_i=\max(0,k_a\delta_i+c_a\dot\delta_i)
\]

転がり抵抗:

\[
F_{rr,i}=-C_{rr}N_i\tanh(v_{x,i}/v_\epsilon)
\]

canonical engineering baseline:

- `k_a = 300000 N/m`
- `c_a = 18000 N s/m`
- `Crr = 0.08`
- `v_epsilon = 0.5 m/s`

位置clampは行わず、nose/main/tail contactはunilateral reactionの結果として決まる。

## 6. FEM cable ground contact

FEM free mass point j の cable centerline に対して半径 `d/2` を含む penetration を使う。

\[
\delta_j=\max(0,z_j+d/2-z_g)
\]

node数依存を弱めるため tributary unstrained length `Ltrib` で ground law をscaleする。

\[
N_j=L_{trib,j}\max(0,k'_c\delta_j+c'_c\dot\delta_j)
\]

滑り摩擦:

\[
F_{f,j}=-\mu_cN_j\tanh(v_{x,j}/v_\epsilon)
\]

canonical engineering baseline:

- `k'_c = 100 N/m^2`
- `c'_c = 2 N s/m^2`
- `mu_c = 0.20`
- `v_epsilon = 0.5 m/s`

位置/速度の強制投影は行わない。

## 7. Pre-launch initialization

### 7.1 Aircraft

2500 N nominal cable forceを考慮し、wheel spring reactionからCG-zとpitchを2変数の静的平衡として解く。x方向は物理時間の発進加速度として残すため、initializerではx位置のみ固定する。

### 7.2 Secant

material natural lengthは **1100 m**。initial geometryとequivalent stiffnessから2500 Nとなる `initial_reel_in_m` を派生量として求める。

v0.13 の1000 m baselineの等価軸剛性を長さ比で維持するengineering conversionとして、dynamic Secant canonicalは

\[
k=20000\frac{1000}{1100}=18181.818\ \mathrm{N/m}
\]

\[
c=2500\frac{1000}{1100}=2272.727\ \mathrm{N\,s/m}
\]

を使用する。この `k,c` はFEMの物性 `EA` と同一視しない。

### 7.3 FEM

- material natural length: `1100 m`
- canonical node count: `20`
- hook/winch fixed during initializer
- cable gravity ON
- cable ground contact ON
- cable aero OFF
- free-node velocity = 0

静的ポテンシャル最小化により cable shape を求め、winch-side tensionが2500 Nとなるよう initial reel-inを反復修正する。

initializer metadata:

- initial CG height / pitch
- initial reel-in
- initial hook/winch tension
- maximum free-node static residual
- initial maximum cable-ground penetration

## 8. State vector

base ODE state は v0.14 で **10 → 9 states** とする。

```text
x_m
z_m
u_m_s
w_m_s
theta_rad
q_rad_s
drum_angle_rad
drum_rate_rad_s
winch_integral
```

全canonical caseで `pilot.full_climb.eas_ki=0` であったため、将来用に残っていた `pilot_integral`、pilot anti-windup、`eas_ki` は削除する。

FEMではこのbase stateに可変個数のfree-node position/velocityを連結する。

## 9. Winch model

Dynamic winch canonical ground-start state:

- drum angle = 0
- drum rate = 0
- initial reel-in はinitializer-derived quantity

FEM final reference element:

\[
L_{N,0}=L_{ref}(N_{removed}+1)-L_{reel,0}-R_D\theta_D
\]

Secant effective natural length:

\[
L_0=L_{unstrained}-L_{reel,0}-R_D\theta_D
\]

full-throttle drum-shaft torque curveはactual club winch powertrain未同定のengineering placeholderであり、1200 ftへfitさせるためには変更しない。

## 10. Standard winch tension control

Prescribed tension および Standard dynamic tension-control target は simulation start ではなく LIFTOFF を基準にする。

- GROUND_ROLL: Clubはlaunch-drive rampを継続。Standard tension-controlは従来どおりpretension targetから開始
- LIFTOFF後: initial rise rate
- PROGRESSIVE_ROTATION後: rotation rise rate
- upper launch: cable angleでrelease targetへ低減

Standard dynamic winchでは `tension_control_drive_torque()` を用い、available full-throttle torqueを超えない。

## 11. Standard pilot control

Standardは `standard_pilot_control()` が担当し、phaseは次の5段階。

```text
GROUND_ROLL
→ INITIAL_CLIMB
→ PROGRESSIVE_ROTATION
→ FULL_CLIMB
→ RELEASE_PREPARATION
```

- initial climb: pitch-attitude feedback
- rotation start: EAS / EAS-rate / AGL条件
- progressive rotation: half-cosine blend
- full climb: EAS feedback + pitch-rate damping
- release preparation: cable-ground-angleでleveling側へblend
- release: low tension / cable angle / pitch rate / pitch reference

Standardはv0.14でもEAS full-climb controlを維持する。

## 12. Club pilot control

Clubは `club_pilot_control()` が担当し、Standardと共有controllerへ押し込まない。

```text
GROUND_ROLL_TAILWHEEL
    ↓ natural liftoff
POST_LIFTOFF_HOLD_DOWN
    ↓ AGL = 60.96 m (200 ft)
PROGRESSIVE_ROTATION
    ↓ transition complete
FULL_CLIMB
    ↓ terminal winch fade
LEVEL_OFF
    ↓ release criteria
RELEASE
```

### 12.1 Ground roll: tailwheel-load control

目的は固定舵角保持ではなく、索が引かれて尾輪が接地した後に正のtailwheel loadを維持すること。

canonical engineering controller:

- elevator bias: `-6.0 deg`
- tailwheel load cue: `500 N`
- feedback gain: `10 deg`
- pitch-rate gain: `2.4 s`

500 Nは実測尾輪荷重ではなくcontroller targetである。

### 12.2 Liftoff → 200 ft: taut-cable hold-down

特定pitch角は追従しない。

- forward elevator: `+0.9 deg`
- low-tension cue: `4000 N`
- low-tension relief: `20 deg`
- end: **AGL 60.96 m**

張力十分なら頭を押さえ、張力低下時だけ前舵を緩める。

### 12.3 Rotation

Rotation開始条件はClubでは **AGL 60.96 mのみ** とする。EAS start threshold / minimum EAS-rateは使用しない。

canonical:

- target: `45 deg`
- duration: `7.4 s`
- rotation pitch Kp: `1.5`

45°はcommand targetであり、modelが必ず45°実姿勢に達するよう空力・舵角をfitしない。

### 12.4 Full climb

Club full climbの主制御はEASではなくnominal pitch attitude。

canonical common Secant/FEM pilot:

- nominal target: `45 deg`
- pitch Kp: **`0.50`**
- elevator bias: **`-4.0 deg`**
- pitch-rate gain: `2.4 s`

この値はSecant/FEMで別々にfitせず、両plantに同じClub SOPを適用した局所探索で、FEM stall marginをguardrailから離すために採用した。

## 13. Club Vw protection

速度はprimary targetではなくupper constraintとして扱う。

### 13.1 Pilot pitch protection

- start margin: `Vw - 1.5 m/s`
- pitch increment gain: `4 deg/(m/s)`
- max pitch increment: `10 deg`

Vwへ接近すると45°commandにpositive pitch incrementを加える。

### 13.2 Pilot → winch power-reduction request

さらにVwへ接近すると `WINCH_POWER_REDUCTION_REQUEST` をlatchする。

canonical request margin: **`Vw - 0.9 m/s`**。

canonical request marginは0.9 m/sを維持する。v0.17のlaunch-power ramp導入後に `wind_x=-4.5 m/s` を再実行し、request 25.993 s → reduced throttle 26.493 s → normal release 45.843 sを確認した。Vw warningは200 ft未満に1回発生するが、200 ft以降のmax EASは32.693 m/s < 33.333 m/sでstop eventはない。一方minimum modeled stall marginは約0.236 m/sまで低下するため、このcaseはfeedback-path verificationでありoperational wind envelopeのvalidationではない。

requestは連続ODE stateにせず `pilot_state.power_reduction_request_time_s` としてaccepted stepでlatchする。

### 13.3 Operator response

- operator response delay: `0.5 s`
- reduced throttle: `0.70`

\[
t\ge t_{request}+0.5\ \mathrm{s}
\]

でspeed-protection reductionを有効化する。

## 14. Club terminal fade

Vw対策のpower reductionと曳航終盤のterminal fadeは別物である。

`club_winch_throttle()` は、

1. pretension-equilibrium commandからlaunch commandへの3.5 s ramp
2. speed-protection reduction
3. terminal fade

を一つの実関数で扱う。LIFTOFFはこの制御則へ入力しない。

simulation `t=0` はpretension成立済みのALL OUT / launch startとし、初期equivalent commandは

\[
 u_0 = \frac{T_{pre} r_D + \tau_{loss}}{M_{full}(0)}
\]

で与える。canonicalでは `T_pre=2500 N`, `r_D=0.5 m`, `tau_loss=120 N m`, `M_full(0)=2880 N m` なので `u_0 ≈ 0.476`。その後、

\[
 u(t)=u_0+\min(t/3.5,1)(u_{launch}-u_0)
\]

とする。`launch_throttle=1.0` はactual pedal/lever positionの測定値ではなく、equivalent full-throttle torque curveに対するprovisional launch-drive targetである。

terminal fade canonical:

- hook→winch LOS ground angle `60 deg` で開始
- `75 deg` でcomplete
- final throttle `0`

FEM local cable tangentは曲率・振動の影響が大きいため、terminal fade geometryは **hook→winch line-of-sight (LOS)** を使う。

fade progressはaccepted-step discrete stateで単調latchし、機体がwinchを通過して瞬間LOS角が戻ってもprogressを逆戻りさせない。

## 15. LEVEL OFF / Club release

terminal fade開始後、Club pilotは45° pitch controlから **ground flight-path angle** controlへ滑らかに移る。

\[
\gamma_{ground,target}=0^\circ
\]

canonical:

- `level_off.target_gamma_ground_deg = 0`
- `gamma_kp = 1.35`
- elevator bias = `0 deg`

Club normal release criteria:

- pilot phase = `LEVEL_OFF`
- terminal fade complete
- cable tension `<= 500 N`
- hook→winch LOS angle `>= 77.5 deg`
- `|q| <= 1.5 deg/s`
- `|gamma_ground| <= 5 deg`

旧Clubの `pitch_reference_deg=8.4°` は使用しない。

## 16. Vw event policy

Clubについてphase-name依存のearly-Vw policyを使用しない。

### AGL < 60.96 m

- Vw crossing: warning
- `1.1 Vw`: hard analysis stop

### AGL >= 60.96 m

- Vw: hard constraint

200 ftをまたぐ1 accepted step内もlinear interpolationし、200 ft crossing時点ですでにEAS>Vwなら、そのcrossing時刻で `MAX_WINCH_SPEED_EXCEEDED` stopとする。

Standardは従来のearly phase policyを維持する。

## 17. LIFTOFF event

最後のaircraft wheel contactのsigned depthが正から負へcrossingした時刻をlinear interpolationして `LIFTOFF` eventとする。

LIFTOFF後にのみ:

- stall stop
- aerodynamic alpha validity stop
- rotation start judgement
- post-liftoff tension/throttle behavior

を有効化する。

## 18. FEM reeling/remesh

RK4 step中はstate dimensionを固定し、accepted step後にwinch-side spanが

\[
l_n<kL_{ref}
\]

となった場合のみ最後の自由質点を削除する。

remesh前後で:

- remaining node positions/velocities
- free reference length
- free cable mass

を保存する。

`node_removal_length_factor=0.2` はengineering baseline。

## 19. 数値積分

- fixed-step RK4
- Prescribed/Secant canonical `dt = 0.01 s`
- FEM canonical `dt = 0.005 s`
- event crossing: candidate step内linear interpolation
- remesh: accepted step endのみ

FEMについて `dt / dt/2 / dt/4`、20/22/40-node sensitivityをvalidation対象とする。

## 20. Diagnostics

共通:

- aircraft ground contact count
- max tire penetration / signed depth
- total normal force / rolling resistance
- nose/main/tail normal load
- cable tension / cable ground angle / LOS angle
- EAS / alpha / stall margin
- gamma_ground / pitch rate
- pilot phase / pitch target
- LIFTOFF time

Club追加:

- `pilot_speed_protection_pitch_increment_deg`
- `pilot_power_reduction_requested`
- `club_speed_reduction_active`
- `winch_terminal_fade_progress`

FEM追加:

- contact node count
- max cable-ground penetration
- ground normal/friction force
- friction power
- free reference length/mass
- active/removed node count
- `cable_node_2...N_x_m / z_m`：可視化・CSV出力用の各free-node位置。winch側で削除済みのnodeは `NaN`

可視化APIは `make_winch_launch_animation(history, fps=..., mode=...) -> fig, ani, axes` とする。Prescribed/Secantはhook→winchを1線分、FEMはhook→free nodes→winchを各物理spanごとの線分として描画し、解析で得たたるみ形状を人工的な補間なしで表示する。動画時間はシミュレーション時間に一致する1×実時間とし、`fps` は再生倍率ではなく時間分解能として扱う。frame数は原則 `round(duration * fps)` とし、MP4/GIF保存時も生成時と同じ `fps` を用いる。x-z表示はfollow / overviewとも物理縮尺1:1（x 1 m = z 1 m）とする。表示範囲は発航を読むために必要な物理領域だけを指定し、overviewの通常の地面側下限は -20 m とする。縮尺1:1は表示範囲の手計算ではなく Matplotlib の equal aspect によって最終レンダリング時に直接保証する。

## 21. Case hierarchy

```text
cases/
├─ g103a_prescribed_tension.json
├─ g103a_dynamic_winch.json
├─ g103a_club_sop_dynamic_winch.json
└─ fem/
   ├─ g103a_dynamic_winch_fem.json
   └─ g103a_club_sop_dynamic_winch_fem.json
```

5ケースを同一1100 m setupで比較し、SOP差とlaunch-distance差を混在させない。

## 22. Canonical representative controls

### 22.1 Prescribed tension / Standard

- ground-roll elevator: `-5.2 deg`
- pitch-rate gain: `2.85 s`
- initial target / bias / Kp: `20 deg / -0.9 deg / 0.80`
- rotation: `24.5 m/s / 6.0 s / 40 deg / Kp 1.1`
- full-climb target EAS / Kp: `28.7 m/s / 0.60`
- tension rise: `0 / 90 N/s`
- tension cap: `6500 N`
- release fade: `58 → 70 deg`

### 22.2 Dynamic Secant / Standard

- pitch-rate gain: `2.40 s`
- initial: `20 deg / -1.5 deg / Kp 0.84`
- rotation: `26.5 m/s / 7.4 s / 40 deg / Kp 1.5`
- full-climb: `27.8 m/s / EAS Kp 0.70`
- target-tension rotation rise: `560 N/s`
- target-tension cap: `6200 N`
- cable: `L0=1100 m`, `k=18181.818 N/m`, `c=2272.727 N s/m`

### 22.3 Dynamic FEM / Standard

- pitch-rate gain: `2.40 s`
- initial: `20 deg / -1.5 deg / Kp 0.84`
- rotation: `26.5 m/s / 7.4 s / 40 deg / Kp 1.5`
- full-climb: `27.5 m/s / EAS Kp 0.75`
- target-tension rotation rise: `560 N/s`
- target-tension cap: **`5000 N`**
- FEM: `1100 m`, `20 nodes`

### 22.4 Club common Secant/FEM SOP

- pitch-rate gain: `2.4 s`
- tailwheel cue / gain: `500 N / 10 deg`
- hold-down: `+0.9 deg`, cue `4000 N`, relief `20 deg`, until `60.96 m`
- rotation: `45 deg`, `7.4 s`, Kp `1.5`
- full climb: `45 deg`, **Kp `0.50`, bias `-4.0 deg`**
- Vw pitch-protection margin/gain: `1.5 m/s / 4 deg/(m/s)`
- winch request margin: `0.9 m/s`
- launch throttle: `1.0` (provisional)
- launch-power ramp duration: `3.5 s`
- operator delay / reduced throttle: `0.5 s / 0.70`
- terminal fade LOS: `60 → 75 deg`

## 23. Representative selection guardrails

- normal release
- no stop event
- no runtime validation error
- airborne minimum modeled stall margin `>= 1.5 m/s`
- maximum absolute pitch rate `<= 8 deg/s`
- weak link `7000 N`未満
- full-climb Vw未超過
- Club Rotation starts at AGL 60.96 m
- no practical pre-rotation cable slack

Clubの実運用「約1200 ft」は重要なvalidation datumだが、optimizer targetとしてtorque curveを高度へfitしない。

## 24. V&V requirements

### Verification

- wheel unilateral contact direction
- FEM cable ground-contact direction / tributary scaling
- pre-launch equilibrium residual / pretension
- FEM removal conservation
- full-throttle drive torque bound
- LIFTOFF event
- Standard/Club controller separation
- Club 200 ft transition
- Club LEVEL OFF + LOS release
- Club Vw early-warning / late-hard-stop
- exact 200 ft crossing boundary
- pilot request → operator delay → power reduction
- normal release

### Numerical V&V

要求:

- full-launch FEM `dt=.005/.0025/.00125`
- 20/22/40-node comparison
- mesh comparison時は、time-step errorとmesh errorを混同しない共通の十分細かいdtを使用する
- half/double contact stiffness with damping scaled by `sqrt(k)`
- cable-ground friction sensitivity

v1.0ではStandard / Club FEMのfull-launch 3段階time-step convergenceを完了した。20-node canonicalでの結果は次のとおり。

| case | dt=.005 | dt=.0025 | dt=.00125 |
|---|---:|---:|---:|
| FEM Standard release altitude | 322.4141 m | 322.4218 m | 322.4257 m |
| FEM Standard max tension | 5.148449 kN | 5.148443 kN | 5.148443 kN |
| FEM Club release altitude | 354.4857 m | 354.4813 m | 354.4764 m |
| FEM Club max tension | 6.301138 kN | 6.301144 kN | 6.301140 kN |

全runでnormal release、stop eventなし。canonical `dt=.005 s` はfull-launch quantitiesについて十分に収束している。

mesh comparisonは20/22/40 nodesすべて共通 `dt=.0025 s` で実施した。

| case | nodes | release altitude | max tension |
|---|---:|---:|---:|
| FEM Standard | 20 | 322.4218 m | 5.148443 kN |
| FEM Standard | 22 | 322.2001 m | 5.148949 kN |
| FEM Standard | 40 | 321.6750 m | 5.149168 kN |
| FEM Club | 20 | 354.4813 m | 6.301144 kN |
| FEM Club | 22 | 354.4528 m | 6.300885 kN |
| FEM Club | 40 | 354.3843 m | 6.298289 kN |

20→40 nodesの差は、Standardでrelease altitude `-0.232%` / max tension `+0.014%`、Clubでrelease altitude `-0.027%` / max tension `-0.045%`。40-nodeはさらに `dt=.00125 s` でも再計算し、両SOPとも `.0025 s` と実質一致した。

したがってv0.17で `40 nodes + dt=.005 s` から得られた大きなpeak-tension差は、時間離散化と空間離散化が混在した結果であり、mesh-only evidenceとして使用しない。v1.0のmesh conclusionは共通 `dt=.0025 s` の結果に基づく。

accepted-step node removal前後のwinch tension jumpは、converged runsではsingle-digit N程度である。remaining node positions/velocities、free reference length、free cable massを保存する既存remesh ruleは維持し、追加のsmoothing state/controllerは導入しない。

開発途中のFEM ground-contact vectorizationで、penetrationが0の節点にもvelocity dampingが作用する差分を検出した。v1.0では旧scalar loopと同じく `penetration > 0` の節点だけcontact dampingを評価し、このunilateral semanticsをunit testで固定する。

v1.0 numerical records:

- `validation/reference/canonical_regression.json`
- `validation/reference/fem_convergence.json`
- `validation/run_fem_convergence.py`（full matrix raw recordを必要時に再生成）
- v0.17 launch-power transient recordsはhistorical/change-specific referenceとして保持

### Physical validation pending

- tire loaded radius / load-deflection
- grass rolling resistance
- actual cable-ground friction
- actual Club ALL OUT→liftoff throttle/RPM/drum-speed history
- actual winch torque-speed / power curve
- aerodynamic ground effect

v0.14 validation files are retained only as historical references for the pre-v0.17 launch-power transient.

## 25. Canonical limitations

- ground contact is compliant penalty contact, not rigid complementarity;
- tire/grass/cable-ground parameters are engineering baselines;
- equivalent winch torque curve is not yet identified from actual club hardware;
- VSPAERO aerodynamic model is local linear and ground effect is not included;
- G72/G73 automatic-release installation-axis mapping is unresolved;
- Club controller loads/cues are engineering representations of SOP, not measurements of pilot neuromuscular response;
- release-local FEM quantities are not claimed fully mesh-converged solely because global altitude/load metrics are stable.

## 26. Release acceptance

v1.0 candidate is considered fully release-accepted only after:

1. all 5 canonical case JSONs load and run;
2. `pytest -q` passes;
3. Notebook code cells compile against the v1.0 package schema;
4. case provenance contains no stale ground_roll_tension_hold / nominal_throttle statements;
5. validation JSON records the final canonical results;
6. clean ZIP extraction repeats `pytest -q` successfully;
7. ZIP integrity test succeeds.
8. FEM animation history contains fixed node-position columns and the animation reduces the segment count after node removal.


## 27. v2.11 change note

Club winch controlの主要なbreaking changeは次の2点である。

- `winch.control.ground_roll_tension_hold` を削除。地上滑走中の2500 N PI holdは行わない。
- `winch.control.nominal_throttle` を `launch_throttle` へ改名し、`launch_power_ramp_duration_s=3.5` を追加。

2500 Nは `launch_setup.pretension_n` とpre-launch equilibriumにのみ残る。launch-power rampはtime-basedで、ground/airborne stateとは独立である。

## 28. v1.0 change note

v1.0ではcase schemaやcanonical physicsを変更していない。主変更は次の3点である。

- FEM force evaluationを配列演算へ整理し、RK4中間段ではhistory diagnosticsを生成しない。
- cable-ground contactのunilateral conditionを `penetration > 0` として明示し、vectorizationで意味が変わらない回帰testを追加する。
- Standard / Club FEMのfull-launch dt/2・dt/4と、共通 `dt=.0025 s` による20/22/40-node mesh V&Vを完了する。

remesh rule、FEM constitutive law、Club launch-power ramp、pilot/winch controllerには新しい補正stateを追加しない。
