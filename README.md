# GaussianPlume — Methane Dispersion Model (C++)

A small C++17 tool that estimates **methane (CH₄) concentration downwind of a
known leak source** using the classic **Gaussian plume model**. Given the GPS
position of a stack/leak, its emission rate, and a set of field sensor readings
(position, measured ppm, wind speed and direction), it computes what the plume
model predicts at each sensor location for every Pasquill–Gifford atmospheric
stability class (A–F), so predictions can be compared against measurements.

The code also contains an (currently disabled) **inverse / reverse model** that
back-calculates the emission rate *Q* from a measured concentration.

> 📐 For the equations, dispersion coefficients and coordinate conventions see
> **[docs/MODEL.md](docs/MODEL.md)**.

---

## Contents

- [What it does](#what-it-does)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Output](#output)
- [Project layout](#project-layout)
- [How the code flows](#how-the-code-flows)
- [Results with the shipped data set](#results-with-the-shipped-data-set)
- [Known issues and limitations](#known-issues-and-limitations)
- [Ideas for future work](#ideas-for-future-work)
- [Credits](#credits)

---

## What it does

```
 Input/inputParams          Input/data.csv
 (source location,          (sensor lat/lon/height, ppm,
  emission rate)             wind speed, wind direction)
          │                          │
          └────────────┬─────────────┘
                       ▼
     lat/lon ──► local x/y metres relative to the source
                       ▼
     rotate each point into the wind-aligned frame
                       ▼
     Gaussian plume C(x, y, z) for stability classes A..F
                       ▼
        Output/forward_model_data.csv
  (measured vs. modelled concentration, g/m³)
```

**Forward model** (active): source + emission rate ➜ predicted concentration at each sensor.

**Reverse model** (implemented in `realm::compute_amountQ`, call site commented
out in `realm::run`): measured concentration ➜ implied emission rate *Q*.

---

## Quick start

### Requirements

- CMake ≥ 3.10
- A C++17 compiler (GCC, Clang or MSVC)
- POSIX threads (linked via `Threads::Threads`)

No third-party libraries are needed.

### Build

```bash
git clone <this-repo> GaussianPlume
cd GaussianPlume
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Run

`GasFlow` takes one optional argument: the project directory that contains
`Input/` and `Output/`. It defaults to `..`, so running from `build/` needs no
argument.

```bash
cd build && ./GasFlow          # from build/ (project dir defaults to ..)
./build/GasFlow .              # or from the repo root
./build/GasFlow /path/to/case  # or any folder with Input/ and Output/
```

The result is written to `<project>/Output/forward_model_data.csv` (overwriting
the existing file). Progress and per-point debug output are printed to stdout;
with the shipped data set the run takes well under a second.

To summarise the result per stability class (Python 3, standard library only):

```bash
python3 scripts/summarize_output.py
```

---

## Inputs

### `Input/inputParams` — source definition

Simple `key = value` file; anything after `//` is a comment.

```text
stack_outlet = (45.529786, 10.174233, 0.42) // (lat, lon, height)
molar_mass   = 16.04                        // (g/mol)
flow_rate    = 56.42                        // (g/hr)
```

| Key            | Meaning                                                         | Units              |
|----------------|-----------------------------------------------------------------|--------------------|
| `stack_outlet` | Leak/stack position as `(latitude, longitude, height)`          | decimal °, °, m    |
| `molar_mass`   | Molar mass of the gas (CH₄ = 16.04). Used to convert sensor ppm to g/m³. **Required, > 0.** | g/mol |
| `flow_rate`    | Emission rate *Q*. Converted internally to g/s (÷ 3600). **Required, > 0.** | g/h |

### `Input/data.csv` — sensor measurements

One row per measurement; the first line is a header and is skipped.

```csv
Latitude,Longitude,Height,Concentration,Wind Velocity,Wind Angle
45.529724,10.174458,0,10.3,0.2,7
```

| Column          | Meaning                                         | Units                 |
|-----------------|-------------------------------------------------|-----------------------|
| `Latitude`      | Sensor latitude (WGS-84)                        | decimal degrees       |
| `Longitude`     | Sensor longitude                                | decimal degrees       |
| `Height`        | Sensor height above ground                      | m                     |
| `Concentration` | Measured CH₄ concentration                      | ppm                   |
| `Wind Velocity` | Wind speed at the time of the measurement       | m/s                   |
| `Wind Angle`    | Direction the wind blows **towards**, counter-clockwise from east (converted to radians on read) | degrees — see [the convention note](#known-issues-and-limitations) |

The parser is positional (column order matters), expects exactly these six
comma-separated columns, and does not tolerate quoted fields or blank lines.

---

## Output

### `Output/forward_model_data.csv`

Each input row produces **six** output rows, one per stability class.

```csv
Latitude,Longitude,Height,Stability,Sensor Measurement,Forward-Model-Concentration
45.529724,10.174458,0.000000,0,0.007416,0.000000
```

| Column                        | Meaning                                                              | Units |
|-------------------------------|----------------------------------------------------------------------|-------|
| `Latitude`, `Longitude`, `Height` | Sensor position, echoed from the input                           | °, °, m |
| `Stability`                   | Pasquill–Gifford class: `0`=A (very unstable) … `5`=F (stable)       | –     |
| `Sensor Measurement`          | Measured value converted from ppm: `ppm × molar_mass / 24.45 / 1000` (ideal gas, 25 °C, 1 atm; ≈ 0.656 mg/m³ per ppm for CH₄) | g/m³  |
| `Forward-Model-Concentration` | Gaussian plume prediction at the sensor. `0` means the sensor is upwind of the source **or** so far off the plume axis that the value underflows. | g/m³ |

To pick the stability class that best fits reality, compare the last two
columns per class (e.g. RMSE across all sensors).

### `Output/data.csv`

Legacy output from an earlier run of the **reverse** model
(`Latitude,Longitude,Height,Stability,Flow-Rate`). It is not regenerated by the
current code and many of its values are numerically meaningless (they blow up
for sensors far off the plume centreline). Kept for reference only.

---

## Project layout

```
GaussianPlume/
├── CMakeLists.txt         # Build definition — compiles src/*.cpp into `GasFlow`
├── Input/
│   ├── inputParams        # Source location, molar mass, emission rate
│   └── data.csv           # 328 field measurements (lat, lon, h, ppm, wind)
├── Output/
│   ├── forward_model_data.csv   # Result of the forward model (generated)
│   └── data.csv                 # Legacy reverse-model output
├── docs/
│   └── MODEL.md           # Physics, equations, coefficients, conventions
├── scripts/
│   └── summarize_output.py  # Per-class summary of the forward-model output
└── src/                   # ← the code that is actually built
    ├── main.cpp           # Input parsing, chunked CSV reading, output writing
    ├── realm.h/.cpp       # Gaussian plume model (ufvm::realm)
    ├── latlon2local.*     # MATLAB-Coder-generated geodetic → ENU conversion (unused)
    ├── cosd.* / sind.*    # Helpers for latlon2local (unused)
    └── CMakeLists.txt     # Copy of the top-level CMakeLists (not used by the build)
```

> **Note:** the `.cpp`/`.h` files at the repository root are an **older copy**
> of `src/` and are **not compiled** (`CMakeLists.txt` only globs `src/*.cpp`).
> The root `realm.h` still has an outdated `compute_Concentration` signature.
> Edit files under `src/` only.

---

## How the code flows

| Step | Function | File |
|------|----------|------|
| 1. Parse source parameters | `readInputParams()` → `setStackOutlet()` | `src/main.cpp` |
| 2. Split `data.csv` into 100-line chunks (byte offsets) | `calculate_offsets()` | `src/main.cpp` |
| 3. For each chunk: parse rows, convert ppm → g/m³ and wind angle ° → rad, build a `ufvm::realm` | `process_chunk()` | `src/main.cpp` |
| 4. Convert each sensor lat/lon to metres east/north of the source (sensor height kept as z) | `realm::distance()` (equirectangular approximation), `realm::run()` | `src/realm.cpp` |
| 5. Rotate into wind-aligned frame | `realm::rotatePointInv()` | `src/realm.cpp` |
| 6. Dispersion widths σᵧ, σ_z for the class & downwind distance | `realm::computeSigmaYSigmaZ()` | `src/realm.cpp` |
| 7. Plume concentration with ground reflection | `realm::compute_Concentration()` | `src/realm.cpp` |
| 8. Concatenate results and write CSV | `main()` | `src/main.cpp` |

Chunks are currently processed **sequentially**. A `std::async` pipeline (up to
5 concurrent chunks) is present but commented out in `main()`; re-enabling it
is safe because each chunk owns its own `realm` instance.

---

## Results with the shipped data set

Run of `Input/data.csv` (328 readings, source at 45.529786 N, 10.174233 E,
0.42 m, Q = 56.42 g/h), summarised with `scripts/summarize_output.py`:

| Class | Non-zero predictions | Max predicted (g/m³) | Median predicted / measured (non-zero rows) |
|-------|---------------------:|---------------------:|--------------------------------------------:|
| A     | 23 / 328 | 4.2 × 10⁻³ | 0.0047 |
| B     | 23 / 328 | 3.6 × 10⁻³ | 0.00013 |
| C     |  9 / 328 | 9.3 × 10⁻⁴ | 0.0038 |
| D     |  3 / 328 | 2.1 × 10⁻⁵ | 7 × 10⁻⁶ |
| E     |  0 / 328 | 0 | – |
| F     |  0 / 328 | 0 | – |

Mean measured concentration: 3.03 g/m³ (median 45.8 ppm, max ≈ 201 000 ppm).

**How to read this**

- The model **under-predicts by 2–6 orders of magnitude**. The best (least bad)
  fit is class A/B, which is consistent with the light winds (0.2–1.4 m/s).
- The sensors are only **2–22 m from the source** (median 4.9 m). The
  Pasquill–Gifford σ curves are calibrated for ≥ 100 m, so at a few metres σᵧ
  and σ_z are centimetres to ~2 m. Off the narrow plume axis the predicted value
  underflows to 0, and in stable classes (E, F) every point does.
- 114 of 328 readings are downwind under the current direction convention, but
  only 23 get a non-zero value for the reason above.
- Very high readings (≫ 10 000 ppm) mean some sensors sat practically at the
  leak — a regime a far-field Gaussian plume cannot represent.

In short, the code now runs consistently end to end, but this data set is
outside the model's validity range. A near-field model or measurements at
≥ 50–100 m would be needed for a meaningful comparison.

## Known issues and limitations

Fixed in the latest revision (no change to the model formulas):

- ✅ Wind angle is now converted from degrees to radians when read.
  Previously degree values were fed straight into `cos`/`sin`.
- ✅ ppm → g/m³ uses `molar_mass / 24.45` instead of a hard-coded 0.72
  (0.72 was the 0 °C value; measured values are now ~9 % lower).
- ✅ Sensor height from the CSV is now passed to the model (previously forced
  to 0; the shipped data is all at 0 m, so results are unaffected).
- ✅ Input/output folder can be passed on the command line.
- ✅ Missing or non-positive `molar_mass` / `flow_rate` now stops the run with
  an error instead of silently producing zeros.

Still open:

1. **Wind direction convention is unconfirmed.** The code treats the angle as
   the direction the wind blows *towards*, counter-clockwise from east. Most
   anemometers report the direction the wind comes *from*, clockwise from
   north; if so, use `θ = 270° − angle`. Checking against the data does not
   settle it: under either convention the median ppm upwind and downwind is
   about the same, which is expected this close to the source in near-calm wind.
2. **Near-field use.** See *Results* — the data is far inside the range where
   the σ curves are valid.
3. **Stability class D coefficient table typo** (`computeSigmaYSigmaZ`): the
   branch after `1000 ≤ dx < 3000` is written as `dx >= 1000 && dx < 10000`
   instead of `dx >= 3000`. Harmless today (the earlier branch wins). Left
   unchanged on purpose — model formulas are not edited.
4. **Class F gaps at bin edges.** Class F uses strict `>` / `<`, so distances
   exactly on a bin edge (200, 700, 1000 m, …) fall back to `a = b = 1`. Also
   left unchanged.
5. **Flat-earth distance.** `realm::distance()` is an equirectangular
   approximation — fine for this data, not for kilometre-scale studies.
6. **Unused / leftover code.** `latlon2local`, `cosd`, `sind`,
   `compute_Concentration_old`, `computeAlpha_x`, several unused `realm`
   members, the stale root-level source copies and the commented block at the
   top of `realm::run()` are remnants of the CFD framework this code came from.
7. **No tests.** Results have not been validated against a reference
   implementation.

## Ideas for future work

- Confirm the anemometer's wind-direction convention.
- Add a regression test against a hand-computed case.
- Re-enable the reverse model and emit `Output/reverse_model_data.csv`;
  combine sensors (e.g. least squares over all points) instead of a per-point
  *Q* estimate.
- Automatically select the best-fit stability class (or estimate it from wind
  speed / solar radiation via the Pasquill table).
- A small Python/Matplotlib script to plot measured vs. modelled values on a map.

