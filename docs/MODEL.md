# Model Reference

This document describes the physics and numerics implemented in
`src/realm.cpp`. It is meant for anyone who wants to check, extend or fix the
model. For build/run instructions see the [README](../README.md).

---

## 1. Gaussian plume equation

For a continuous point source at effective height *H* in a steady, uniform
wind *u* blowing along the +x axis, the concentration at (x, y, z) is

```
              Q                  ⎛   y²   ⎞ ⎡    ⎛ (z − H)² ⎞       ⎛ (z + H)² ⎞ ⎤
C(x,y,z) = ─────────── · exp ⎜ − ────── ⎟·⎢exp ⎜ − ──────── ⎟ + exp ⎜ − ──────── ⎟ ⎥
           2π u σᵧ σ_z       ⎝   2σᵧ²  ⎠ ⎣    ⎝   2σ_z²  ⎠       ⎝   2σ_z²  ⎠ ⎦
```

| Symbol | Meaning | Units | Source in code |
|--------|---------|-------|----------------|
| `Q` | Emission rate | g/s | `flow_rate` from `inputParams` (g/h ÷ 3600) |
| `u` | Wind speed | m/s | `Wind Velocity` column |
| `σᵧ, σ_z` | Lateral / vertical dispersion widths, functions of x and stability | m | `computeSigmaYSigmaZ()` |
| `H` | Source height (no plume rise modelled) | m | 3rd value of `stack_outlet` |
| `x, y, z` | Downwind, crosswind, vertical distance from the source | m | after rotation (§3) |
| `C` | Concentration | g/m³ | `Forward-Model-Concentration` |

The second exponential in the bracket is the **image source** that models total
reflection of the plume at the ground.

Concentration is defined only for `x > 0` (downwind of the source); for
`x ≤ 0` the code returns `0`.

Implemented in `realm::compute_Concentration()`.

### Inverse (reverse) model

Solving the same equation for *Q* gives

```
Q = C_meas · 2π u σᵧ σ_z / [ exp(−y²/2σᵧ²) · ( exp(−(z−H)²/2σ_z²) + exp(−(z+H)²/2σ_z²) ) ]
```

Implemented in `realm::compute_amountQ()`, with optional multiplicative
uncertainty factors `(1 + GPS_unc)` on x/y/z and `(1 + Sensor_unc)` on C. The
call is commented out in `realm::run()`. Note that this per-sensor inversion
becomes ill-conditioned far from the plume centreline (the denominator → 0),
which is why the legacy `Output/data.csv` contains astronomically large values.

---

## 2. Dispersion coefficients (σᵧ, σ_z)

The code uses the **Pasquill–Gifford curves in the ISC3 (EPA Industrial Source
Complex) rural parameterisation**, with downwind distance `x` in **kilometres**:

```
σ_z = a · x^b                    (capped at 5000 m)
θ   = 0.017453293 · (c − d · ln x)      [radians]
σᵧ  = 465.11628 · x · tan θ
```

`c`, `d` are fixed per class; `a`, `b` depend on the class **and** on the
downwind distance band.

| Class (code) | P-G class | Atmosphere | c | d |
|---|---|---|---|---|
| 0 | A | Very unstable | 24.1670 | 2.5334 |
| 1 | B | Moderately unstable | 18.3330 | 1.8096 |
| 2 | C | Slightly unstable | 12.5000 | 1.0857 |
| 3 | D | Neutral | 8.3330 | 0.72382 |
| 4 | E | Slightly stable | 6.2500 | 0.54287 |
| 5 | F | Moderately stable | 4.1667 | 0.36191 |

The `a`, `b` distance bands are coded literally in `computeSigmaYSigmaZ()` and
follow the ISC3 user's guide (Vol. II, Table 1-1). Any class/distance not
covered falls back to `a = b = 1` (and `c = d = 1` for an unknown class).

### Choosing a class

The model runs all six classes for every sensor. As a rough guide (Pasquill,
daytime):

| Surface wind (m/s) | Strong sun | Moderate sun | Slight sun | Night, cloudy | Night, clear |
|---|---|---|---|---|---|
| < 2 | A | A–B | B | E | F |
| 2–3 | A–B | B | C | E | F |
| 3–5 | B | B–C | C | D | E |
| 5–6 | C | C–D | D | D | D |
| > 6 | C | D | D | D | D |

With the very light winds in the shipped data set (0.2–1.4 m/s), classes A–B
(daytime) or E–F (night) are the physically plausible ones.

---

## 3. Coordinate frames

### Geographic → local metres

`realm::distance(lat₁, lon₁, lat₂, lon₂)` uses an equirectangular projection
around the source (R = 6 371 000 m):

```
x_east  = R · (λ₂ − λ₁) · cos((φ₁ + φ₂)/2)
y_north = R · (φ₂ − φ₁)
```

Accurate to well under 1 % for distances of a few kilometres. The source then
sits at `(0, 0, H)`; the sensor's z is its `Height` from the CSV.

### Local → wind-aligned frame

`realm::rotatePointInv(p, source, θ)` rotates a point by −θ about the source:

```
x' =  x·cos θ + y·sin θ      (downwind)
y' = −x·sin θ + y·cos θ      (crosswind)
```

This is correct if θ is the direction the wind blows **towards**, measured
**counter-clockwise from east, in radians**.

The `Wind Angle` column is read in **degrees** and converted to radians in
`process_chunk()` (`src/main.cpp`), so the input is expected to be the
direction the wind blows **towards, counter-clockwise from east**.

> ⚠️ **Convention not yet confirmed.** Most anemometers report the
> meteorological direction (where the wind comes **from**, clockwise from
> north). If that is what the sensor logged, convert it before use:
>
> ```
> θ_math = (270 − wind_from_deg) · π / 180
> ```

---

## 4. Unit conversions

| Quantity | Conversion in code | Comment |
|---|---|---|
| Emission rate | `g/h ÷ 3600 → g/s` | `readInputParams()` |
| Measured concentration | `ppm × M / V_m / 1000 → g/m³` | `process_chunk()`. `M` = `molar_mass` from `inputParams`, `V_m` = 24.45 L/mol (ideal gas, 25 °C, 1 atm). For CH₄: 0.656 mg/m³ per ppm. |
| Wind direction | `° × π / 180 → rad` | `process_chunk()` |

---

## 5. Assumptions

Standard Gaussian plume assumptions apply; results are only meaningful when
they roughly hold:

- Steady-state emission and meteorology over the averaging time.
- Wind speed and direction uniform with height and in space.
- Flat, unobstructed terrain; no buildings or canopy.
- No plume rise (buoyancy/momentum), no deposition, no chemical decay.
- Downwind transport dominates diffusion along x (unreliable for `u ≲ 1 m/s`,
  which applies to much of the shipped data).
- σ curves were derived for ~10-minute averages and distances ≥ 100 m; very
  near-source predictions (tens of metres, as here) are extrapolations.

## References

- U.S. EPA (1995). *User's Guide for the Industrial Source Complex (ISC3)
  Dispersion Models, Vol. II – Description of Model Algorithms*,
  EPA-454/B-95-003b.
- Turner, D. B. (1994). *Workbook of Atmospheric Dispersion Estimates*, 2nd ed.
- Stockie, J. M. (2011). The mathematics of atmospheric dispersion modeling.
  *SIAM Review*, 53(2), 349–372.
