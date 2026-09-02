# winch_launch_simulator implementation specification v1.3.0

## 1. Scope

v1.3.0 retains the v1.2.0 two-axis architecture and changes two physical baselines: pitch inertia and pre-launch support initialization.

1. **plant fidelity**
   - prescribed tension + ideal Secant
   - dynamic winch + distributed FEM cable
2. **pilot / operation policy**
   - RWTH Reference
   - Club

The old lumped elastic `Dynamic Secant` model and the old `Standard SOP` are removed rather than retained as compatibility paths.

## 2. Canonical cases

```text
RWTH Reference / Prescribed Tension + Ideal Secant
Club Pilot     / Prescribed Tension + Ideal Secant
RWTH Reference / Dynamic Winch + FEM
Club SOP       / Dynamic Winch + FEM
```

Case files:

```text
cases/g103a_rwth_reference_prescribed_tension.json
cases/g103a_club_prescribed_tension.json
cases/g103a_rwth_reference_dynamic_winch_fem.json
cases/g103a_club_dynamic_winch_fem.json
```

No additional canonical JSON cases are permitted without a specific validation purpose.

## 3. Coordinate and aircraft state

The basic state remains

```text
[x_G, z_G, u, w, theta, q, theta_D, omega_D, I_W]
```

for dynamic-winch cases. Prescribed cases keep the same base state layout but drum states are inert.

Positive conventions remain those of v1.1.0:

- ground `x`: toward winch
- ground `z`: downward
- positive pitch: nose-up
- positive elevator: trailing-edge down / nose-down aerodynamic moment

## 4. Aerodynamics

`CL` and `Cm` use the VSPAERO stability-derivative model around AoA 2°:

```text
CL = CL_ref + CL_alpha*(alpha-alpha_ref) + CL_qhat*qhat + CL_deltae*deltae
Cm = Cm_ref + Cm_alpha*(alpha-alpha_ref) + Cm_qhat*qhat + Cm_deltae*deltae
```

Drag uses

```text
CD = CD0 + K*CL^2
```

with

```text
CD0 = 0.007906026873656432
K   = 0.025021597731889444
```

from the `CLwtot/CDwtot` static polar fit.

## 5. Ideal Secant plant

The ideal Secant contains no physical cable state.

```text
mass        = 0
elasticity  = 0
damping     = 0
gravity     = 0
drag        = 0
ground load = 0
```

Cable force direction is the hook-to-winch line of sight. Its magnitude is supplied by `prescribed_tension()`.

The old fields

```text
stiffness_n_m
damping_n_s_m
unstrained_length_m  # for a secant spring
```

must not appear in ideal-Secant cases.

## 6. Common prescribed-tension excitation

Both prescribed cases use exactly the same open-loop input. It is an **analysis excitation**, not either operational SOP.

LIFTOFF-relative points:

```text
0.0 s   2500 N
1.5 s   1600 N
18.0 s  4750 N
```

Linear interpolation is used between points. Terminal fade is based on hook-to-winch LOS ground angle:

```text
60°: no terminal reduction yet
75°: 400 N
```

The prescribed pair therefore supports a controlled pilot-policy comparison under identical external tension forcing.

## 7. FEM cable plant

Canonical FEM:

```text
natural length             1100 m
node count                 20
diameter                   0.005 m
linear density             0.0188 kg/m
Young's modulus            4.025e10 Pa
strain-rate stress coeff.  1.56e7 Pa s
cable CD                   1.3
weight                     enabled
aerodynamic drag           enabled
ground interaction         enabled
reel-in                    enabled
```

The FEM equations, WLO07 cylindrical drag completion, ground model, and node-removal mechanism are unchanged from v1.1.0.

## 8. Dynamic winch hardware

Both FEM cases use the same plant:

```text
drum radius          0.5 m
equivalent inertia   20 kg m²
loss torque           120 N m
```

and the same equivalent full-throttle drum-shaft torque curve.

This curve remains an engineering baseline until the actual Club winch powertrain is identified.

## 9. RWTH Reference pilot

Source-derived structure:

- elevator trim `-3°`
- trim-only below `50 m`
- target EAS `30 m/s`
- airspeed/dynamic-pressure feedback after safety altitude
- pitch-rate damping in parallel
- second-order fade-in after crossing safety altitude

Implementation:

```text
speed_kp          1.5
pitch_rate_gain   2.4 s
fade_tau1         1.0 s
fade_tau2         2.0 s
```

The source papers describe the controller form but do not provide the final human-controller gains. These numerical values are therefore G103A analysis settings, not literature constants.

## 10. RWTH dynamic-winch policy

The dynamic FEM RWTH case uses tension control.

Target structure is source-derived:

- moderate initial force
- controlled rise after liftoff
- maximum target force
- terminal reduction

G103A-adapted numerical settings:

```text
pretension        2500 N
rise              500 N/s after LIFTOFF
cap               5000 N
terminal fade     cable ground angle 65° → 75°
final target       400 N
```

PI actuator settings:

```text
Kp  0.5 N m / N
Ki  0.02 N m / (N s)
```

Drive torque is clipped by the same full-throttle torque-speed capability used by the Club plant.

## 11. Club pilot

The same `pilot` block must be present in both Club cases.

### Ground roll

```text
elevator bias        -6°
tailwheel cue         500 N
tailwheel gain         10°
```

### Hold-down

```text
until AGL             60.96 m
forward elevator      +0.9°
low-tension cue       4000 N
relief gain             20°
```

The tension cue is a surrogate controller input, not a measured direct pilot tension indication.

### Progressive rotation

```text
duration              7.4 s
nominal target          45°
pitch Kp               1.5
elevator bias           -5°
```

### Full climb

Canonical v1.3.0 common settings:

```text
nominal target          45°
pitch Kp               0.25
elevator bias           -2°
pitch-rate gain         2.4 s
Vw pitch protection    starts at Vw - 1.5 m/s
power-down request     at Vw - 0.9 m/s
```

These gains are analysis settings selected so that the *same Club controller* remains viable on both plant fidelities. They are not copied from the deleted Standard controller.

### Level off and release

Terminal power/tension fade drives the transition into level-off. Level-off controls

```text
gamma_ground → 0°
```

rather than directly imposing a zero pitch attitude.

Normal Club release requires terminal fade completion, low tension, sufficient LOS angle, and the configured pitch-rate / flight-path-angle guardrails.

## 12. Club dynamic-winch policy

```text
launch throttle                   1.0
launch-power ramp                 3.5 s
operator delay after Vw request   0.5 s
reduced throttle                  0.70
terminal fade LOS                 60° → 75°
final throttle                    0
```

The launch-power ramp is meaningful only for Dynamic winch + FEM. It does not exist in the prescribed-tension pair.

## 13. Operational limits and analysis thresholds

Canonical model values remain:

```text
mass                 580 kg
Vw                   33.3333 m/s model-EAS threshold
Vs,1g                20.8333 m/s model-EAS mapping
weak-link threshold  7000 N analysis value
hook load limit      11700 N
```

The Flight Manual operational speeds are not claimed to be source-labelled EAS values; the simulator maps them to EAS because pitot/static/instrument errors are not separately modeled.

## 13A. Pitch inertia baseline

The canonical pitch inertia is

```text
Iyy = 2000 kg m²
```

This is a derived analysis value, not a measured aircraft mass-property result. AMC 22.481 permits `iy = 0.225 LR` when `iy` cannot be determined by more rational means, with `LR` the fuselage overall length without rudder. Using the Maintenance Manual fuselage length `8.18 m` as the `LR` proxy and `m=580 kg` gives `iy=1.8405 m` and `Iyy=1964.7 kg m²`; the simulator rounds this to 2000 kg m².

## 14. Initialization

All four cases start from the same ground geometry and 2500 N pretension condition.

Ideal Secant:

- no reel-in state is needed
- initial reel-in metadata = 0

FEM:

- static cable shape is solved
- initial reel-in is adjusted until winch-side pretension is 2500 N

The aircraft portion of the initializer uses a mainwheel+tailwheel two-point static support state under 2500 N pretension. Nosewheel contact is forbidden and checked explicitly. With the canonical geometry and compliant tire model:

```text
rigid two-point pitch         5.49924 deg
compliant equilibrium pitch  5.25148 deg
CG height                     0.580987 m
mainwheel reaction            5683.52 N
tailwheel reaction               4.944 N
nosewheel reaction               0 N
nosewheel clearance              0.130828 m
```

The FEM cable is then relaxed from that aircraft pose. Representative FEM initialization is approximately:

```text
initial reel-in        4.03253 m
winch tension          2499.94 N
hook tension           2499.97 N
```

## 15. Numerical integration

Fixed-step classical RK4.

```text
Ideal Secant  dt = 0.01 s
FEM           dt = 0.005 s
```

Discrete topology changes occur only after an accepted RK4 step.

### v1.3.0 numerical checks

Club FEM early-transient, canonical 3.5 s power ramp:

| dt [s] | LIFTOFF [s] | distance [m] | EAS [m/s] |
|---:|---:|---:|---:|
| .005 | 4.770926 | 59.12986 | 26.48703 |
| .0025 | 4.770968 | 59.13096 | 26.48720 |
| .00125 | 4.771001 | 59.13187 | 26.48732 |

RWTH FEM full-launch spot check:

```text
dt=.005   release altitude 241.3443 m, max tension 5306.475 N
dt=.0025  release altitude 241.3377 m, max tension 5306.501 N
```

The full v1.1.0 20/22/40-node mesh matrix remains relevant to the unchanged FEM plant equations, but v1.3.0 does not relabel those old trajectories as current-policy results.

## 16. Canonical regression v1.3.0

| case | release altitude [m] | max EAS [m/s] | max tension [kN] | min stall margin [m/s] | max |q| [deg/s] |
|---|---:|---:|---:|---:|---:|
| RWTH / Prescribed | 347.581 | 31.335 | 4.750 | 3.926 | 6.023 |
| Club / Prescribed | 345.334 | 33.058 | 4.750 | 1.988 | 5.359 |
| RWTH / FEM | 241.344 | 32.560 | 5.306 | 3.666 | 7.214 |
| Club / FEM | 329.102 | 32.905 | 5.108 | 1.704 | 4.805 |

All four cases reach normal release with no stop event.

## 17. Comparison semantics

Do not claim that every cell differs by only one physical factor.

- Prescribed RWTH vs Prescribed Club: same forcing, different pilot policy
- RWTH Prescribed vs RWTH FEM: same RWTH pilot policy, but prescribed forcing is replaced by dynamic tension-control plant
- Club Prescribed vs Club FEM: same Club pilot policy, but prescribed input is not the real Club throttle law; this is reduced-order comparison only
- Club FEM: primary Club validation case

## 18. Tests and reference records

Current release acceptance:

```text
python -m pytest -q
63 passed
```

Reference files:

```text
validation/reference/canonical_regression.json
validation/reference/launch_power_ramp.json
validation/reference/vw_feedback.json
```

The old v1.1.0 convergence record must not be presented as a v1.3.0 trajectory record.

## 19. Provenance classes

Every nontrivial parameter should be classifiable as one of:

```text
source_confirmed
source_derived_structure
model_mapping
analysis_setting
engineering_baseline
unresolved
```

In particular, do not label controller tuning as source-confirmed merely because the controller architecture comes from a source.

## 20. Removed interfaces

v1.3.0 retains the removals introduced in v1.2.0:

```text
pilot.sop = "standard"
cable.model = elastic/dynamic secant
secant stiffness_n_m
secant damping_n_s_m
standard_pilot_control()
```

No compatibility translation is provided. Old case JSONs should fail validation rather than silently map to the new model.
