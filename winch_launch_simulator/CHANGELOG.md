# Changelog

## v1.0

- GitLab repository layoutへ整理。
- Python実装を `winch_launch_simulator/` packageへ移動。
- canonical caseを `cases/` 直下へ統合。
- `tests/` と数値V&Vの `validation/` を分離したまま維持。
- 実フライトIGCのraw dataを `data/flight_validation/raw/` へ移動。
- Notebookを `notebooks/` へ集約。
- 再生成可能な `outputs/`、`fig/`、IGC派生CSV/PNGをGit管理対象から除外。
- versioned release-note filesを本ファイルへ統合。
- 計算モデル・canonical parameter・数値アルゴリズム自体は変更していない。

---

# winch_launch_simulator v0.18 release notes

v0.18 は v0.17 の Club launch-power ramp と canonical physics を維持しつつ、FEM 計算経路を整理し、未完だった full-launch numerical V&V を完了する release である。

## Implementation changes

- WLO07 cable aerodynamic drag を element-by-element Python loop から NumPy array operation へ整理した。
- FEM internal axial-force accumulation と cable-ground contact evaluation を配列演算へ整理した。
- fixed-step RK4 の k1–k4 中間評価では、history 用 diagnostics dictionary を生成せず derivative のみ計算するようにした。
- accepted point では従来どおり full diagnostics を生成する。
- Club launch-power ramp、pilot SOP、winch controller、FEM constitutive law、reeling/remesh rule は変更していない。

## Vectorization regression fix

開発途中の contact vectorization で、penetration が 0 の節点にも downward velocity damping が作用する差分を検出した。旧 scalar loop は `penetration > 0` の節点だけ contact damping を評価するため、v0.18 では同じ unilateral semantics を明示的な mask で維持する。

この差分を検出する unit test を追加した。

## Full-launch FEM time-step convergence

20-node canonical FEM を `dt=.005/.0025/.00125 s` で再計算した。

| SOP | metric | .005 s | .0025 s | .00125 s |
|---|---|---:|---:|---:|
| Standard | release altitude [m] | 322.4141 | 322.4218 | 322.4257 |
| Standard | max tension [kN] | 5.148449 | 5.148443 | 5.148443 |
| Club | release altitude [m] | 354.4857 | 354.4813 | 354.4764 |
| Club | max tension [kN] | 6.301138 | 6.301144 | 6.301140 |

両 SOP とも normal release、stop event なし。canonical `dt=.005 s` は full-launch quantities について十分に収束している。

## FEM mesh convergence

時間離散化と空間離散化を混同しないため、20/22/40-node の mesh comparison は共通 `dt=.0025 s` で実施した。

| SOP | nodes | release altitude [m] | max tension [kN] |
|---|---:|---:|---:|
| Standard | 20 | 322.4218 | 5.148443 |
| Standard | 22 | 322.2001 | 5.148949 |
| Standard | 40 | 321.6750 | 5.149168 |
| Club | 20 | 354.4813 | 6.301144 |
| Club | 22 | 354.4528 | 6.300885 |
| Club | 40 | 354.3843 | 6.298289 |

20→40 node の差は、Standard で release altitude `-0.232%` / max tension `+0.014%`、Club で release altitude `-0.027%` / max tension `-0.045%`。

40-node はさらに `dt=.00125 s` でも再計算し、Standard / Club とも `.0025 s` と実質一致した。したがって v0.17 の `40 nodes + dt=.005 s` で見えた peak-load anomaly は mesh-only evidence としては使用しない。

## Remesh continuity

converged runs では accepted-step node removal 前後の winch tension jump は single-digit N 程度であり、既存の remaining-node state / free reference length / free cable mass preservation tests も通過している。したがって v0.18 では remesh algorithm 自体へ補正 state や smoothing controller を追加しない。

## Regression

- `pytest -q`: **56 passed**
- 5 canonical cases: **all normal release, warning/stopなし**
- Standard FEM canonical: release altitude **322.414 m**
- Club FEM canonical: release altitude **354.486 m**
- Club FEM canonical max tension: **6.301 kN < 7.0 kN weak-link analysis threshold**

Numerical records:

- `validation/canonical_regression_v0.18.json`
- `validation/fem_convergence_v0.18.json`
- `validation/run_v018_fem_convergence.py` はfull convergence matrixの再生成用runner

## API compatibility

v0.18 で case schema / public simulation API の新しい breaking change はない。v0.17 の breaking changes (`ground_roll_tension_hold` removal, `nominal_throttle -> launch_throttle`) をそのまま継承する。

## Release acceptance

- current tree: `pytest -q` **56 passed**
- 5 canonical regression: **all normal release, warning/stopなし**
- all packaged JSON files parse successfully
- all Python modules compile successfully
- 3 notebooks / 27 code cells compile successfully
- clean ZIP extraction: `pytest -q` **56 passed**
- ZIP integrity test: **pass**

---

# winch_launch_simulator v0.17 release notes

v0.17はClub SOPの発航開始時winch power applicationを修正するreleaseである。

## Breaking changes

- `winch.control.ground_roll_tension_hold` 削除
- `winch.control.nominal_throttle` → `launch_throttle`
- `winch.control.launch_power_ramp_duration_s` 追加（canonical 3.5 s）
- Clubの `target_tension_n` はground roll中もcontroller targetを持たないためNaN
- historyに `winch_launch_power_ramp_progress` 追加

## Control change

旧実装はLIFTOFFまで2500 N pretensionをPI holdし、LIFTOFF後にlaunch throttleへ切り替えていた。v0.17ではsimulation t=0をpretension成立済みのALL OUT / launch startと定義し、pretension平衡drive commandからlaunch commandまで3.5 s linear rampを行う。rampはLIFTOFFをまたいで継続し、LIFTOFFをwinch-control switchとして使用しない。

初期commandは `u0=(Tpre*rD+tau_loss)/M_full(0)` から算出する。canonical値では `u0≈0.476`。

`launch_throttle=1.0` はactual Club winchの物理lever位置を同定した値ではなく、equivalent torque-curve model上のprovisional launch targetである。

既存のVw speed protection、0.5 s operator delay、terminal fade、Club pilot SOP、Secant/FEM plantは維持する。

## V&V

Release candidateで以下を確認した。

- `pytest -q`: **55 passed**
- 5 canonical launch regression: **全ケース normal release、warning/stopなし**
- Secant Club, 3.5 s ramp: liftoff **4.92714 s**, ground distance **64.281 m**
- FEM Club, canonical `dt=0.005 s`: release altitude **354.486 m**
- FEM Club launch transient `dt/dt2/dt4 = 0.005/0.0025/0.00125 s`:
  liftoff **4.761879 / 4.761986 / 4.762000 s**、ground distance **58.9171 / 58.9200 / 58.9203 m**
- Notebook code cells compile against the v0.17 schema

結果は `validation/canonical_regression_v0.17.json` と `validation/launch_power_ramp_vv_v0.17.json` に保存する。FEM full-launch dt/2・dt/4 とmesh sensitivityは、今回変更した初期launch transientの収束確認とは分けてhistorical/追加V&Vとして扱う。

---

# winch_launch_simulator v0.15 release notes

v0.15は、v0.14で確立したClub SOP・1100 m launch setup・FEM cableの数値モデルを変更せず、アニメーションを解析結果と同じ実時間・物理縮尺で読めるように仕上げるreleaseである。

主な変更:

- `make_winch_launch_animation(history, ...) -> fig, ani, axes` へAPIを統一
- reel-in / FEM node removalに応じて表示segment数を減少
- animationを1×実時間化し、`fps`を実時間1秒あたりの表示frame数として定義
- follow / overviewのx-z物理縮尺をMatplotlib `equal` aspectで1:1化
- overviewの通常の地面側下限を -20 m に維持し、不要な地下領域を排除
- 保存MP4/GIFについて実時間長とpixel/m等価性を回帰テスト

v0.14から継承する数値モデル:

- `standard_pilot_control()` / `club_pilot_control()` へ分離
- Club: tailwheel hold → taut-cable hold-down → AGL 200 ft rotation → nominal 45 deg full climb → LEVEL OFF → RELEASE
- Vw proximity pitch protectionと `WINCH_POWER_REDUCTION_REQUEST`
- operator delay後のwinch speed-protection reductionとterminal fadeを分離
- LEVEL OFFを `gamma_ground -> 0` で定義
- Club release / terminal fadeにhook-to-winch LOS angleを使用
- base continuous stateを10→9 statesへ簡素化 (`pilot_integral`削除)
- all canonical launch distance = 1100 m
- Secant equivalent stiffness/dampingを1100 mへ換算
- FEM natural length = 1100 m, canonical node count = 20
- Club Secant/FEM共通pilot candidate: full-climb Kp 0.50, elevator bias -4.0 deg

Final canonical release altitude:

- Prescribed: 362.229 m
- Secant Standard: 368.724 m
- Secant Club: 363.274 m
- FEM Standard: 322.414 m
- FEM Club: 338.476 m

Numerical V&V:

- FEM Standard dt=.005/.0025/.00125: 322.4141 / 322.4218 / 322.4257 m
- FEM Club dt=.005/.0025/.00125: 338.4764 / 338.4667 / 338.4647 m
- 20→40-node release-altitude change at dt=.0025: Standard -0.232%, Club -0.031%
- final contact-stiffness and cable-ground-friction sensitivities are stored in `validation/club_sop_vv_v0.14.json`

Validation records:

- `validation/canonical_regression_v0.14.json`
- `validation/club_sop_vv_v0.14.json`

Known physical-validation gaps remain the actual winch powertrain, tire/grass parameters, actual cable-ground friction, aerodynamic ground effect, and automatic-release installation-axis mapping.

Numerical note:

- v0.15 does not change the flight-dynamics, controller, winch, Secant, or FEM numerical model from v0.14.
- Therefore `validation/canonical_regression_v0.14.json` and `validation/club_sop_vv_v0.14.json` remain the canonical numerical evidence for v0.15.
