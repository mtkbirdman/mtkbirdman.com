# winch_launch_simulator v1.3.0

G103A の縦運動ウインチ発航を、地上静止から LIFTOFF、曳航上昇、level off、release まで計算する研究用シミュレーターです。

v1.3.0 は v1.2.0 の4ケース architecture を維持し、G103A pitch inertia と pre-launch ground initialization の provenance を改善します。`Iyy` は根拠のない 2500 kg m² から、AMC 22.481 の `iy=0.225 LR`、`LR=8.18 m`、`m=580 kg` による 1964.7 kg m² を丸めた **2000 kg m²** へ変更しました。索張り合わせ時の機体は、主輪＋尾輪の2点接地静的平衡から開始し、前輪が非接触であることを検証します。

## Canonical architecture

| Plant | RWTH Reference | Club |
|---|---|---|
| Prescribed tension + Ideal Secant | `g103a_rwth_reference_prescribed_tension.json` | `g103a_club_prescribed_tension.json` |
| Dynamic winch + FEM cable | `g103a_rwth_reference_dynamic_winch_fem.json` | `g103a_club_dynamic_winch_fem.json` |

### 1. Prescribed tension + Ideal Secant

低忠実度 plant です。

- cable mass = 0
- cable elasticity / damping = 0
- cable gravity = 0
- cable aerodynamic drag = 0
- cable-ground interaction = 0
- hook と winch を直線で結ぶ
- prescribed tension の大きさを、その直線方向へ作用させる

Gäb & Santel の Secant model と同じ理想化を採用します。旧 v1.1.0 の `k_C`, `c_C` を持つ弾性 Secant は完全に削除しました。

Prescribed pair の2ケースは **同じ張力履歴**を使用します。したがって低忠実度 plant 上での差は pilot policy の差です。

canonical open-loop excitation は LIFTOFF 基準で、

```text
0.0 s   2500 N
1.5 s   1600 N
18.0 s  4750 N
terminal fade: LOS 60° → 75°, 400 N
```

です。この履歴は RWTH SOP でも Club SOP でもなく、比較用の共通解析入力です。

### 2. Dynamic winch + FEM cable

高忠実度 plant です。

- dynamic drum inertia
- equivalent full-throttle torque-speed curve
- distributed cable mass
- axial elasticity and damping
- gravity
- cable aerodynamic drag
- cable-ground contact / friction
- reel-in and node removal

を連成します。

RWTH case は tension-control、Club case は throttle-based operation を使用します。FEM cable と winch hardware は2ケースで共通です。

## Pilot / operation policies

### RWTH Reference

Santel (2008) / Gäb & Santel (2011) の reference procedure から、次を source-derived structure として実装します。

- safety altitude 50 m までは trim-only
- trim elevator = -3°
- safety altitude 後に airspeed / dynamic-pressure control と pitch damper を fade-in
- target EAS = 30 m/s
- winch force は初期に抑え、LIFTOFF 後に増加させ、上限を設け、終盤に低減

公開資料には human-controller gain の具体値がないため、G103A 用 controller gain と fade constants は **analysis settings** です。

Dynamic FEM の target tension は、G103A へ適合させた解析設定として、

```text
initial       2500 N
rise rate      500 N/s after LIFTOFF
cap           5000 N
terminal fade  LOS 65° → 75°, 400 N
```

を用います。ASK 21 reference case の数値をそのまま G103A に移植したものではありません。

### Club

Club policy は plant fidelity に依存せず同一の pilot model を使います。

```text
GROUND_ROLL_TAILWHEEL
  ↓ natural LIFTOFF
POST_LIFTOFF_HOLD_DOWN to 60.96 m AGL
  ↓
PROGRESSIVE_ROTATION to nominal 45°
  ↓
FULL_CLIMB nominal 45° + Vw protection
  ↓
LEVEL_OFF: gamma_ground → 0°
  ↓
RELEASE
```

canonical full-climb settings は、両 plant で共通に

```text
pitch target       45°
elevator bias      -2°
pitch Kp           0.25
pitch-rate gain    2.4 s
```

です。これらの数値は current G103A model に対する analysis settings であり、訓練資料の規定値ではありません。

Dynamic FEM では Club winch を throttle operation としてモデル化します。

```text
launch throttle ramp     3.5 s → 1.0
Vw protection request    Vw - 0.9 m/s
operator delay           0.5 s
reduced throttle         0.70
terminal fade            LOS 60° → 75° → throttle 0
```

## Aircraft aerodynamics

`CL` と `Cm` は `G103A.stab` の local-linear stability derivatives を使用します。

抗力は v1.1.0 以降、

```text
CD = CD0 + K * CL^2
```

です。

`G103A.polar` の `CLwtot/CDwtot`, AoA -2°..12° を最小二乗 fit した値は、

```text
CD0 = 0.007906026873656432
K   = 0.025021597731889444
R²  = 0.9999565364889669
```

です。

Repository には再現用として、

- `aero_reference/G103A.polar`（VSPAERO original static polar）
- `aero_reference/G103A_drag_polar.csv`
- `aero_reference/G103A_drag_polar_fit.json`
- `validation/fit_drag_polar.py`

を含みます。

## Pre-launch two-point support

4ケースは 2500 N pretension が成立した **主輪＋尾輪の2点接地静的平衡**から開始します。初期化では主輪・尾輪だけを active contact として垂直力・ピッチモーメント平衡を解き、前輪が地面に接触していないことを明示的に検証します。

canonical geometry / tire model では、

```text
rigid two-point pitch         5.4992 deg
compliant equilibrium pitch  5.2515 deg
CG height                     0.58099 m
mainwheel reaction            5683.5 N
tailwheel reaction               4.94 N
nosewheel reaction               0 N
nosewheel clearance             0.1308 m
mainwheel penetration           18.95 mm
tailwheel penetration            0.0165 mm
```

となります。尾輪反力は小さいものの正であり、2点接地条件を満たします。タイヤを完全剛体として接触深さゼロから開始すると compliant ground model の反力がゼロになるため、time-domain ODE へ連続に接続する目的で静的タイヤ圧縮を含む平衡姿勢を使用します。

## Pitch inertia

AMC 22.481 は、より合理的な方法で pitch radius of gyration を決定できない場合、

```text
iy = 0.225 LR
```

を用いることを認め、`LR` を rudder を除く fuselage overall length としています。G103A の解析値 `LR=8.18 m` と canonical mass `580 kg` を使うと、

```text
iy  = 1.8405 m
Iyy = m iy² = 1964.7 kg m²
```

となるため、canonical model では **`Iyy = 2000 kg m²`** に丸めます。これは実機の直接実測値ではなく、AMC-based derived analysis value です。

## Numerical integration

時間積分は fixed-step classical RK4 です。

- Ideal Secant canonical: `dt = 0.01 s`
- FEM canonical: `dt = 0.005 s`

FEM node removal は accepted RK4 step の終了後だけ実施し、RK4 substage 中の state dimension は固定します。

v1.2.0 では architecture/control change 後に Club dynamic-FEM の最初の6 sについて、

```text
dt = .005 / .0025 / .00125 s
```

を再計算しています。canonical 3.5 s throttle ramp の liftoff time は、

```text
4.770926 / 4.770968 / 4.771001 s
```

で、liftoff distance と EAS も実質一致します。

また RWTH dynamic-FEM では full-launch `dt=.005` と `.0025` の release altitude が、

```text
241.3443 m / 241.3377 m
```

で一致することを確認しました。

旧 v1.1.0 の 20/22/40-node full-launch mesh study は Git history / archived reference として保持します。v1.3.0 では FEM constitutive / discretization equations は変更していませんが、Iyy と初期姿勢を変更したため current-policy time-step checks は再生成しています。

## v1.3.0 canonical regression

| case | release altitude [m] | max EAS [m/s] | max tension [kN] | min stall margin [m/s] | max |q| [deg/s] |
|---|---:|---:|---:|---:|---:|
| RWTH / Prescribed + Ideal Secant | 347.581 | 31.335 | 4.750 | 3.926 | 6.023 |
| Club / Prescribed + Ideal Secant | 345.334 | 33.058 | 4.750 | 1.988 | 5.359 |
| RWTH / Dynamic winch + FEM | 241.344 | 32.560 | 5.306 | 3.666 | 7.214 |
| Club / Dynamic winch + FEM | 329.102 | 32.905 | 5.108 | 1.704 | 4.805 |

4ケースとも normal release、stop event なしです。

この4ケースを **同一SOPの単純な cable-model sensitivity** と解釈してはいけません。特に dynamic pair では winch policy も RWTH force-control / Club throttle と異なります。比較目的は次のように分けます。

- Prescribed pair: same forcing 下での pilot-policy comparison
- RWTH prescribed → RWTH FEM: reduced plant と dynamic distributed plant の method comparison
- Club prescribed → Club FEM: reduced-order comparison。prescribed case は実 Club winch operation の再現ではない
- Club FEM: 実 Club trajectory validation の主対象

## Running

```bash
python run_simulation.py cases/g103a_rwth_reference_prescribed_tension.json
python run_simulation.py cases/g103a_club_prescribed_tension.json
python run_simulation.py cases/g103a_rwth_reference_dynamic_winch_fem.json
python run_simulation.py cases/g103a_club_dynamic_winch_fem.json
```

Regression:

```bash
PYTHONPATH=. python validation/run_regression.py
```

Club early transient / launch-power-ramp:

```bash
PYTHONPATH=. python validation/run_launch_power_ramp.py
```

Vw feedback path:

```bash
PYTHONPATH=. python validation/run_vw_feedback.py
```

Tests:

```bash
python -m pytest -q
```

v1.3.0 release candidate: **63 passed**.

## Provenance status

Source-confirmed / source-derived:

- G103A Flight Manual operational limits used in the model
- G103A geometry items explicitly traced to dimensioned manual drawings
- VSPAERO `.stab`, `.vspaero`, `.vsp3`, static polar data
- RWTH reference pilot structure and force-control philosophy
- Santel/Gäb FEM cable concept and material properties used for the reference synthetic cable

Analysis settings / unresolved hardware-specific values include:

- `Iyy = 2000 kg m²`: AMC 22.481 の `iy=0.225 LR` と `LR=8.18 m`, `m=580 kg` から得る 1964.7 kg m² を丸めた analysis value
- tailwheel vertical coordinate `z = +0.061 m`
- effective loaded mainwheel radius `0.220 m`
- equivalent winch torque-speed curve
- Club controller numerical gains
- RWTH G103A controller numerical gains
- prescribed-tension common excitation
- deterministic weak-link break threshold `7000 N`

The simulator is a research model, not an operational approval or flight-training authority.
