# winch_launch_simulator v1.0

G103A の縦運動ウインチ発航を、**地上静止状態から ground roll・LIFTOFF・曳航上昇・level off・release まで**計算する研究用シミュレーターです。現在の基準実装は v1.0、実装仕様は `docs/implementation_spec.md`、開発履歴は `CHANGELOG.md` です。


## Repository layout

```text
winch_launch_simulator/
├─ README.md
├─ CHANGELOG.md
├─ pyproject.toml
├─ .gitignore
├─ .gitlab-ci.yml
├─ run_simulation.py
├─ winch_launch_simulator/   # solver implementation
│  └─ models/                # physical submodels
├─ cases/                    # canonical inputs
├─ tests/                    # pytest correctness tests
├─ validation/               # numerical V&V runners and references
├─ notebooks/                # simulation/search/flight analysis
├─ data/flight_validation/raw/
└─ docs/
```

`outputs/`、図・動画、IGC派生CSVは再生成可能なのでGit管理しません。Python実装をpackage化したため、Notebook・tests・validation runnerで `sys.path` を書き換える必要はありません。

## Install

```bash
python -m pip install -e .
```

テスト・Notebookまで含める場合は、

```bash
python -m pip install -e ".[dev,notebook]"
```

## 現在のモデル構成

同一の ground-start setup から、次の3段階で物理忠実度を上げます。

1. **Prescribed tension**
   - 機体3自由度 + wheel-ground contact
   - tension historyを外部規定
2. **Dynamic winch + Secant cable**
   - drum dynamics / drive torque / winch control
   - 索を等価ばね・ダンパで表現
   - Standard / Club SOP
3. **Dynamic winch + FEM cable**
   - distributed elastic cable
   - cable mass / weight / aerodynamic drag
   - cable-ground contact / friction
   - reel-in / accepted-step node removal
   - Standard / Club SOP

Standard / Club は同じplantへ異なる操縦・ウインチ操作を適用するcaseとして扱います。

## Canonical case tree

```text
cases/
├─ g103a_prescribed_tension.json
├─ g103a_dynamic_winch.json
├─ g103a_club_sop_dynamic_winch.json
├─ g103a_dynamic_winch_fem.json
└─ g103a_club_sop_dynamic_winch_fem.json
```

通常のcanonical FEM baselineは20 nodes / `dt=0.005 s`です。高解像度確認では40 nodes / `dt=0.0025 s`を使用でき、v1.0のmesh/time-step V&Vで両者のglobal launch quantitiesが十分近いことを確認しています。

## Canonical physical setup

```text
launch distance       1100 m
pre-launch pretension 2500 N
weak-link analysis    7000 N
Vw                     33.3333 m/s (120 km/h)
```

G103A representative geometry, body-fixed +x forward / +z down:

```text
belly hook      (+0.535, +0.366) m
nosewheel       (+2.029, +0.508) m
mainwheel       (-0.169, +0.366) m
tailwheel       (-4.537, +0.061) m
```

Aircraft ground baseline:

```text
normal stiffness   300000 N/m
normal damping      18000 N s/m
rolling Crr             0.08
velocity smoothing      0.5 m/s
```

FEM cable-ground baseline:

```text
normal stiffness       100 N/m²
normal damping           2 N s/m²
sliding friction mu      0.20
velocity smoothing       0.5 m/s
```

これらのground parameterはengineering baselineであり、field validation対象です。

## Pre-launch initializer

`winch_launch_simulator/initialization.py` がphysical `t=0`を作ります。simulation開始時点はpretension成立済みの `ALL OUT / launch start` とします。

```text
measured aircraft geometry
        ↓
wheel vertical/pitch equilibrium
        ↓
1100 m material cable
        ↓
FEM static equilibrium with gravity + ground contact
        ↓
derived initial reel-in
        ↓
winch tension ≈ 2500 N
        ↓
physical simulation t = 0
```

## Club SOP

Club caseでは、ground rollからreleaseまで次のphaseを使用します。

```text
GROUND_ROLL_TAILWHEEL
→ natural LIFTOFF
→ POST_LIFTOFF_HOLD_DOWN to 200 ft
→ PROGRESSIVE_ROTATION toward 45°
→ FULL_CLIMB nominal 45° + Vw protection
→ terminal winch fade
→ LEVEL_OFF (gamma_ground → 0°)
→ RELEASE
```

Secant/FEMで同じpilot/winch commandを使用します。

```text
pitch-rate gain              2.4 s

ground roll:
  elevator bias             -6.0 deg
  tailwheel load cue         500 N
  tailwheel feedback gain     10 deg

hold-down:
  until AGL                  60.96 m
  forward elevator           +0.9 deg
  low-tension cue            4000 N
  low-tension relief           20 deg

rotation:
  target                       45 deg
  duration                    7.4 s
  pitch Kp                    1.5

full climb:
  nominal target               45 deg
  pitch Kp                    0.50
  elevator bias              -4.0 deg
  speed-protection margin     1.5 m/s
  pitch increment gain        4 deg/(m/s)
  power-reduction request     Vw - 0.9 m/s

winch:
  launch throttle             1.00 (provisional)
  launch-power ramp            3.5 s
  operator delay              0.5 s
  reduced throttle            0.70
  terminal fade LOS           60 → 75 deg
```

`launch_throttle=1.0` はequivalent full-throttle torque curveに対するcanonical launch-drive targetです。actual Club winchのlever position、engine/transmission/drum特性は今後のhardware identification対象です。

Club launch commandはpretension平衡に必要なequivalent commandから始まり、3.5 sのlinear rampで`launch_throttle`へ到達します。このrampはground rollとLIFTOFFをまたいで連続し、LIFTOFFをwinch-controlの切替条件には使用しません。

## Secant 1100 m baseline

1000 m equivalent baselineの軸方向特性を長さ比で1100 mへ換算し、

```text
unstrained length  1100 m
stiffness          18181.818 N/m
damping             2272.727 N s/m
```

を使用します。この `k,c` はFEMのmaterial property `EA` と同一視しません。

## Canonical regression

v1.0の最終release treeから5 canonical casesを再実行した結果です。

| case | release altitude [m] | max EAS [m/s] | max tension [kN] | min airborne stall margin [m/s] | max |q| [deg/s] |
|---|---:|---:|---:|---:|---:|
| Prescribed | 362.229 | 33.269 | 5.154 | 2.881 | 7.935 |
| Secant Standard | 368.724 | 32.288 | 6.236 | 2.158 | 7.897 |
| Secant Club | 377.174 | 31.653 | 6.263 | 2.498 | 5.471 |
| FEM Standard | 322.414 | 30.461 | 5.148 | 1.952 | 7.942 |
| FEM Club | 354.486 | 32.656 | 6.301 | 1.563 | 5.542 |

Club ground run:

| case | liftoff [s] | ground distance [m] | liftoff EAS [m/s] | mean ground acceleration [m/s²] |
|---|---:|---:|---:|---:|
| Secant Club | 4.927 | 64.28 | 26.17 | 5.31 |
| FEM Club | 4.762 | 58.92 | 26.49 | 5.56 |

両caseとも `t=0` のequivalent commandは `0.475694...`、3.5 sで `launch_throttle=1.0` へ到達します。詳細時系列は `validation/reference/launch_power_ramp.json` に保存しています。

## FEM numerical V&V

20-node FEMをfull launchで `dt=.005/.0025/.00125 s` まで再計算しました。

| SOP | metric | .005 s | .0025 s | .00125 s |
|---|---|---:|---:|---:|
| Standard | release altitude [m] | 322.4141 | 322.4218 | 322.4257 |
| Standard | max tension [kN] | 5.148449 | 5.148443 | 5.148443 |
| Club | release altitude [m] | 354.4857 | 354.4813 | 354.4764 |
| Club | max tension [kN] | 6.301138 | 6.301144 | 6.301140 |

全runでnormal release、stop eventなしです。canonical `dt=.005 s` のglobal full-launch quantitiesは十分に収束しています。

mesh sensitivityは20/22/40 nodesすべて共通 `dt=.0025 s` で比較しました。

| SOP | nodes | release altitude [m] | max tension [kN] |
|---|---:|---:|---:|
| Standard | 20 | 322.4218 | 5.148443 |
| Standard | 22 | 322.2001 | 5.148949 |
| Standard | 40 | 321.6750 | 5.149168 |
| Club | 20 | 354.4813 | 6.301144 |
| Club | 22 | 354.4528 | 6.300885 |
| Club | 40 | 354.3843 | 6.298289 |

20→40 nodesの差はStandardでrelease altitude `-0.232%` / max tension `+0.014%`、Clubでrelease altitude `-0.027%` / max tension `-0.045%`です。40-nodeは `dt=.00125 s` でも再計算し、`.0025 s` と実質一致しています。

詳細は `validation/reference/fem_convergence.json` を参照してください。

## Vw feedback verification

ClubはAGL 200 ft未満のVw crossingをwarningとして扱い、200 ft以降はVwをhard constraintとします。200 ft crossingをまたぐ1 integration stepも補間します。

`wind_x=-4.5 m/s` のSecant Club regressionでは、

```text
Vw proximity
→ pilot pitch protection
→ WINCH_POWER_REDUCTION_REQUEST  t=25.993 s
→ 0.5 s operator delay
→ reduced throttle active        t=26.493 s
→ normal release                 t=45.843 s
```

を確認しています。200 ft以降のmax EASは `32.693 m/s < Vw 33.333 m/s` です。minimum modeled stall marginは約 `0.236 m/s` まで低下するため、このcaseはfeedback pathのverificationに使用し、validated operational wind envelopeとは扱いません。

## Physical validation status

数値V&Vとは別に、2026-04-05 / 11 / 12 / 19の実フライト12発航のIGCを使ってClub caseの軌跡を比較しました。1 Hz GNSSログのためwheel-offやreleaseを秒以下で同定する用途には使わず、ground-track accelerationと発航軌跡のenvelopeを主に見ています。

- ground rollからGNSS liftoff proxyまでの時間・距離は現在の3.5 s launch-power rampと矛盾しない
- zero-wind caseは実軌跡より浅い上昇となる
- 3～4 m/s程度のheadwindを与えると、一部の実フライトで20～30 sの水平距離・高度が大きく改善する
- 実風、winch RPM、drum speed、cable tensionが未取得のため、pilot controllerやwinch torque curveの再fitには使用していない

実フライト比較はパッケージ内のrelease acceptance recordには含めず、physical validationの継続課題として扱います。

## Run

```bash
python run_simulation.py cases/g103a_prescribed_tension.json
python run_simulation.py cases/g103a_dynamic_winch.json
python run_simulation.py cases/g103a_club_sop_dynamic_winch.json
python run_simulation.py cases/g103a_dynamic_winch_fem.json
python run_simulation.py cases/g103a_club_sop_dynamic_winch_fem.json
```

結果は `outputs/<case_id>/` へ保存されます。

動画を保存する場合は、

```bash
python run_simulation.py cases/g103a_club_sop_dynamic_winch_fem.json --animation
```

とします。`make_winch_launch_animation(history, fps=..., mode=...)` は計算済み `history` を直接受け取り、`fig, ani, axes` を返します。Prescribed / Secantはhook→winchの1線分、FEMは保存された各節点位置を結ぶ実要素線分として描画します。動画は1×実時間で、x-z表示はfollow / overviewとも物理縮尺1:1です。MP4/GIF保存APIは `save_mp4(ani, path=..., fps=..., dpi=...)` / `save_gif(...)` です。

## Notebook

Main notebook:

```text
notebooks/winch_launch_simulator.ipynb
```

Parameter search:

```text
notebooks/winch_launch_parameter_search.ipynb
notebooks/winch_launch_fem_parameter_search.ipynb
```

Club search schemaは `ground_roll / hold_down / rotation / full_climb / level_off / club_throttle` を使用します。

## Tests

```bash
python -m pytest -q
```

v1.0 release treeでは **56 passed** を確認しています。

主要coverage:

- ground-start / LIFTOFF
- Standard / Club controller separation
- tailwheel-load feedback
- Club 200 ft transition
- Club Vw early-warning / late-hard-stop
- exact 200 ft boundary interpolation
- LEVEL OFF + LOS release
- club throttle mode
- power-reduction request / operator delay / activation
- winch torque bound
- FEM constitutive / ground contact / reeling / node-removal conservation
- FEM node-position history / element-by-element animation / release hide
- case schema / provenance

## Validation records

現行branchでは、現在の基準値だけをversionlessなreference名で保持します。過去版の差分はGit historyと `CHANGELOG.md` から追跡します。

```text
validation/
├─ run_regression.py
├─ run_fem_convergence.py
└─ reference/
   ├─ canonical_regression.json
   ├─ fem_convergence.json
   ├─ launch_power_ramp.json
   └─ vw_feedback.json
```

`run_regression.py` は5 canonical casesを再実行して `canonical_regression.json` を生成します。FEMのfull convergence matrixは `run_fem_convergence.py` で再生成します。

v1.0ではsolver numericsを変更していないため、初期repositoryでは実在するv0.18 numerical recordをcurrent referenceへ引き継ぎ、各JSONの `source_record` / `source_version` に出典を残しています。以後、solver physics・numerical scheme・canonical parameterを変更したreleaseではrunnerを再実行してcurrent referenceを更新します。

## Release acceptance

v1.0について次を確認済みです。

- 5 canonical case JSONがload/run可能
- 5 canonical casesがnormal release、warning/stopなし
- `python -m pytest -q`: 56 passed
- Python modules compile
- Notebook code cells compile
- validation JSON parse
- legacy Club winch-control schema referenceなし
- clean ZIP extraction後のpytest通過
- ZIP integrity check通過

## Known limitations

- ground contactはcompliant penalty spring/damperです。
- tire loaded radius / grass rolling resistance / cable-ground frictionは未同定です。
- Club launch-power ramp 3.5 sとtailwheel cueはengineering controller modelです。
- `launch_throttle=1.0` とfull-throttle winch torque curveはactual hardware未同定のparameterです。
- aerodynamic ground effectは含みません。
- VSPAERO aerodynamic modelは局所線形モデルで、大舵角域はmodel validity確認が必要です。
- G72/G73 automatic-release installation-axis mappingは未確定です。
- FEM release-local quantityはglobal altitude/loadが安定していても完全mesh-convergedとは扱いません。
- flight-by-flight physical validationには発航地点の実風、機体実重量、winch RPM/drum speed、cable tensionが必要です。
