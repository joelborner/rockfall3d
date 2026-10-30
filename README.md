# rockfall3d

[![Build and Test](https://github.com/joelborner/rockfall3d/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/joelborner/rockfall3d/actions/workflows/build-and-test.yml)

`rockfall3d` is a C++ rigid-body rockfall simulator for batch trajectory simulations on raster DEMs with roughness.

## Contact
Joël Borner

WSL Institute for Snow and Avalanche Research SLF

joel.borner@slf.ch

## Features

- DEM-based terrain interaction with bilinear interpolation
- rigid-body dynamics with quaternion orientation
- nonsmooth dynamics with hard contact laws, Coulomb friction, and optional per-contact rolling resistance. Moreau timestepping scheme and JOR/SOR-prox contact solve for terrain interaction
- roughness model based on Perlin noise
- deterministic random orientation and roughness, reproducible with seeds
- shared-memory, thread-based parallel execution
- trajectory and/or raster output
- optional terrain export and debugging options

## Repository Layout

- [`rockfall3d/`](rockfall3d): source code
- [`sample_data/`](sample_data): example XML files and sample inputs
- [`CMakeLists.txt`](CMakeLists.txt): CMake build configuration
- [`CMakePresets.json`](CMakePresets.json): CMake presets for the tested Windows build
- [`interface/`](interface): executable-provided parameter catalogue embedded during the build
- [`docs/interface-description.md`](docs/interface-description.md): versioned interface contract for GUI integrations
- [`tests/`](tests): numerical tests and executable interface-command checks
- [`vcpkg.json`](vcpkg.json): vcpkg dependency manifest used by local builds and CI
- [`.github/workflows/`](.github/workflows): GitHub Actions workflows for automated build and test runs
- [`vcpkg-triplets/`](vcpkg-triplets): custom vcpkg triplets used for optimized CI/release builds
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md): bundled third-party attribution
- [`REFERENCES.md`](REFERENCES.md): conceptual references for important methods used in the code

## Requirements

### Bundled in the repository

- Quickhull by Anatoliy V. Tomilov, bundled in [`rockfall3d/external/quickhull/include/quickhull.hpp`](rockfall3d/external/quickhull/include/quickhull.hpp)

### External requirements

- CMake 3.21 or newer
- a C++17 compiler
- Eigen 3.4.x
- GDAL

## Platform

- currently tested on Windows
- Linux and macOS are not officially supported yet

## Installation on Windows

Install [vcpkg](https://github.com/microsoft/vcpkg) and bootstrap it:

```powershell
git clone https://github.com/microsoft/vcpkg.git
.\vcpkg\bootstrap-vcpkg.bat
```

The repository uses `vcpkg.json` manifest mode. When CMake is configured with the vcpkg toolchain file, project dependencies are resolved automatically from the manifest.

The tested build path uses the Microsoft C++ toolchain (MSVC). If Eigen is installed manually instead of through vcpkg, configure CMake with `-DEIGEN3_INCLUDE_DIR=...`.

## Build

Using the included preset:

```powershell
$env:VCPKG_ROOT="C:\path\to\vcpkg"
cmake --preset windows-release
cmake --build --preset windows-release
```

Equivalent manual commands:

```powershell
cmake -S . -B build\windows-release `
  -G "Visual Studio 17 2022" -A x64 `
  -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT\scripts\buildsystems\vcpkg.cmake"

cmake --build build\windows-release --config Release
```

If CMake does not find Eigen automatically, add `-DEIGEN3_INCLUDE_DIR=...` to the configure command.

The executable is written to:

```text
build\windows-release\Release\rockfall3d.exe
```

## Running

To run a simulation, pass the XML config file:

```powershell
.\build\windows-release\Release\rockfall3d.exe path\to\config.xml
```

GUI integrations can query the selected executable's supported parameter sections:

```powershell
.\build\windows-release\Release\rockfall3d.exe --describe-interface
```

This prints an embedded JSON description without running a simulation or creating
files. The initial catalogue covers optional rock models only; this branch does
not support fragmentation. See the [interface contract](docs/interface-description.md)
for partial coverage and compatibility behaviour. Interface-command tests run
through CTest using CMake; the plugin project owns form and XML-mapping tests.

## Input Files Referenced by XML

### Rock list file

The file referenced by `<rocks list="..."/>` contains one rock point-cloud path (`.pts`) per line.

- relative paths are resolved relative to the list file location

### Start points file

The file referenced by `<starts list="..."/>` contains:

```text
X Y Z dz
```

where `dz` is optional:

- if 3 columns are defined, `X`, `Y`, `Z` are the absolute start coordinates of the centre of mass of the rock.
- if 4 columns are defined,`dz` is the desired initial clearance between the terrain and the lowest rock vertex after the initial orientation is applied. In this case, `Z` is ignored. 

## XML Configuration Reference

The parser expects self-closing tags, except that `<terrain>` may be a block
when terrain interactions are configured.

Required tags:

- `<terrain/>`
- `<rocks/>`
- `<starts/>`
- `<simulation/>`
- `<stop/>`
- `<integrator/>`
- `<solver/>`

Optional tags:

- `<roughness/>`
- `<output>...</output>`
- `<debug/>`

The top-level `<rockfall3d/>` marker is accepted in sample files, but the parser does not rely on it.

### `<terrain/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `dem` | Path to the DEM raster (`.tif`, `.tiff`, or `.asc`) |

Without additional interaction laws, the existing self-closing form remains valid:

```xml
<terrain dem="data/dem_2m.tif"/>
```

Terrain interaction laws are optional self-closing tags inside a terrain block:

```xml
<terrain dem="data/dem_2m.tif">
  <interaction type="coulomb_friction" mu="0.5"/>
  <interaction type="rolling_friction" mu_r="0.05"/>
</terrain>
```

Currently supported interaction types:

| Type | Attributes | Meaning |
|---|---|---|
| `coulomb_friction` | `mu` | Constant non-negative Coulomb friction coefficient |
| `coulomb_friction` | `mumin`, `mumax`, `kappa`, `beta` | Slippage-dependent Coulomb friction according to Leine et al. (2014) |
| `rolling_friction` | `mu_r` | Non-negative dimensionless terrain coefficient |

Exactly one Coulomb parameterization is allowed. For slippage-dependent friction,
`mumax` must be at least `mumin`, `kappa` has units `m^-1`, and `beta` has units
`s^-1`. The coefficient is

```text
mu(s) = mumin + (mumax - mumin) * (2/pi) * atan(kappa * s).
```

The per-trajectory state `s` is the accumulated center-of-mass travel distance
during active terrain contact. A contact is active when its solved normal impulse
is positive. Outside active contact, `s` decays exponentially according to
`ds/dt = -beta*s`, so `mu` approaches `mumin` continuously. Coulomb and rolling
friction can be enabled together as independent parts of the same contact law.

The rolling law is the five-dimensional, per-contact formulation of Acary and
Bourrier (2021). As in PlatRock, the solver converts the terrain coefficient to
the dimensional rolling-resistance length `mu_r * equivalent rock radius`, which
bounds the rolling moment impulse. Omitting the interaction, or setting `mu_r="0"`,
selects the original Coulomb-only solver path. The value must be calibrated for
the terrain; `0.05` is illustrative, not a universal default.

### `<rocks/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `list` | Path to the rock list file (`.txt`) |
| `density` | Rock density in `kg/m^3` |
| `random_orientations` | Number of random orientations per rock/start point |

### `<starts/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `list` | Path to the start point list file (`.txt`) |

### `<simulation/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `dt` | Time step |
| `output_dir` | Output directory |

Optional attributes:

| Attribute | Default if omitted | Meaning |
|---|---|---|
| `dump_step` | falls back to `dt` | Output sampling interval |
| `en` | `0` | Normal restitution coefficient, must be in `[0,1]` |
| `et` | `0` | Tangential restitution coefficient, must be in `[0,1]` |
| `parallel_threads` | `1` | Number of parallel worker threads requested |
| `random_seed` | random device | Fallback seed for both orientation and roughness if dedicated seeds are omitted |
| `orientation_seed` | `random_seed` if provided, otherwise random device | Seed for random orientations |
| `roughness_seed` | `random_seed` if provided, otherwise random device | Seed for roughness realization |

For backward compatibility, constant friction may instead be supplied as `mu`
on `<simulation/>` when no `coulomb_friction` terrain interaction is present.

Notes:

- If `parallel_threads="0"` is given, it is clamped to `1`.
- The code clamps the requested thread count to the number of available logical processors.
- If `en` and `et` are omitted, the solver uses the current fully inelastic impact law with zero restitution.

### `<stop/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `max_time` | Hard upper bound on simulated time |
| `max_steps` | Hard upper bound on accepted steps |
| `settle_min_steps` | Number of steps before settle checks become active |
| `settle_linear_speed` | Translational stop threshold in `m/s` |
| `settle_angular_speed` | Angular stop threshold in `rad/s` |

### `<solver/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `type` | Solver type. Current options: `jor-prox`, `sor-prox` |
| `alpha` | Relaxation parameter `0 < alpha < 2` |
| `atol` | Absolute tolerance |
| `rtol` | Relative tolerance |
| `maxiter` | Maximum iterations per contact solve |

### `<integrator/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `type` | Rotational integrator. Options: `schweizer`, `capobianco` |

The `schweizer` integrator uses the energy-consistent pseudo-implicit rotational update described by Schweizer (2015). The `capobianco` integrator denotes the complete published stability-preserving method of Leine et al. (2021): implicit midpoint free rotation followed by the spin-preserving orientation correction. This is the method named `CAPOBIANCO2` in the predecessor code.

Optional attributes for `capobianco`:

| Attribute | Default | Meaning |
|---|---|---|
| `atol` | `1e-10` | Absolute tolerance for the free-rotation Newton solve |
| `rtol` | `1e-10` | Relative tolerance for the free-rotation Newton solve |
| `maxiter` | `50` | Maximum Newton iterations per time step |

References: [Schweizer (2015)](https://doi.org/10.3929/ethz-a-010464319) and [Leine et al. (2021)](https://doi.org/10.1007/s11044-021-09792-y).

### `<roughness/>`

This tag is optional. It is only needed if roughness is enabled.

Required attributes when roughness is enabled:

| Attribute | Meaning |
|---|---|
| `raster` | Roughness raster whose positive finite cell values are roughness grid sizes |

Cells with value `0`, nodata, or NaN disable roughness locally. Query points outside the roughness raster extent are also treated as no roughness.

Optional attributes:

| Attribute | Default if omitted | Meaning |
|---|---|---|
| `step_func` | `smoothstep` | Interpolation method used in the roughness model. Current options: `linear`, `smoothstep`, `smootherstep` |

Notes on roughness raster:
 - Grid size of the roughness raster should be at least 2x bigger than its highest value.
 - Values in roughness raster should be smaller than DEM resolution

### `<output>...</output>`

This tag is optional.

If omitted:

- output falls back to text trajectory files

If present:

- it may contain at most one `<trajectory .../>`
- it may contain one or several `<raster .../>`
- it may contain at most one `<shape .../>`

#### `<trajectory/>`

This child tag is optional.

Required attributes:

| Attribute | Meaning |
|---|---|
| `format` | Output format: `binary` or `text` |
| `metric` | Comma-separated trajectory metrics (see table below) |

#### `<raster/>`

This child tag is optional. Multiple `<raster .../>` tags are allowed.

Required attributes:

| Attribute | Meaning |
|---|---|
| `metric` | Raster metric to aggregate (see table below) |
| `resolution` | Output raster resolution |

Optional attributes:

| Attribute | Default if omitted | Meaning |
|---|---|---|
| `statistic` | required except for `metric="count"` | Raster statistic: `min`, `max`, `sum`, `mean`, `median`, `percentile` |
| `p` | not used | Percentile value when `statistic="percentile"` |
| `name` | auto-generated from metric/statistic | Output raster file stem |

Notes:

- all `<raster .../>` tags inside one `<output>` block must currently use the same `resolution`

#### `<shape/>`

This child tag is optional.

Required attributes:

| Attribute | Meaning |
|---|---|
| `metric` | Shape output type. Current option: `deposition` |

Valid metric names / aliases:

| Variable | Meaning | Trajectory | Raster | Shape |
|---|---|---|---|---|
| `t` | time [s] | X |  |  |
| `x` | X coordinate [m] | X | X |  |
| `y` | Y coordinate [m] | X | X |  |
| `z` | Z coordinate [m] | X | X |  |
| `quat` | Orientation (Quaternion) | X |  |  |
| `vx,vy,vz` | Translational velocities [m/s] | X | X |  |
| `vtot` | Resultant translational velocity [m/s] | X | X |  |
| `wx,wy,wz` | Rotational velocities [rad/s] | X | X |  |
| `wtot` | Resultant rotational velocity [rad/s] | X | X |  |
| `ekin_trans` | Kinetic translational energy [J] | X | X |  |
| `ekin_rot` | Kinetic rotational energy [J] | X | X |  |
| `ekin_tot` | Kinetic energy [J] | X | X |  |
| `jumpH` | Jump height [m] | X | X |  |
| `count` | Number of contributing trajectories per raster cell |  | X |  |
| `deposition` | Final accepted deposition point per trajectory |  |  | X |

### `<debug/>`

This tag is optional.

Optional attributes:

| Attribute | Default if omitted | Meaning |
|---|---|---|
| `terrain_export` | not written | Path for exported terrain raster; must be paired with `terrain_export_resolution` |
| `terrain_export_resolution` | not written | Resolution for exported terrain raster; must be paired with `terrain_export` |
| `timing_report` | `false` | Write `timing_report.txt` to `output_dir` |

Notes:

- `terrain_export` and `terrain_export_resolution` must either both be present or both be absent.

## Example XML

```xml
<rockfall3d/>

<terrain dem="data/dem_2m.tif">
  <interaction type="coulomb_friction"
               mumin="0.2" mumax="1.3" kappa="0.7" beta="20"/>
  <interaction type="rolling_friction" mu_r="0.05"/>
</terrain>

<rocks
    list="data/rocks.txt"
    density="2500"
    random_orientations="100"/>

<starts list="data/start_points.txt"/>

<simulation
    dt="0.002"
    dump_step="0.02"
    en="0.0"
    et="0.0"
    output_dir="output"
    parallel_threads="20"
    orientation_seed="12345"
    roughness_seed="67890"/>

<stop
    max_time="200.0"
    max_steps="1000000"
    settle_min_steps="100"
    settle_linear_speed="0.1"
    settle_angular_speed="0.1"/>

<integrator type="schweizer"/>

<solver
    type="jor-prox"
    alpha="0.5"
    atol="1e-8"
    rtol="1e-6"
    maxiter="10000"/>

<roughness
    raster="data/roughness_grid_size.tif"
    step_func="smoothstep"/>

<output>
  <trajectory
      format="binary"
      metric="t,x,y,z,quat,vx,vy,vz,vtot,wx,wy,wz,wtot,ekin_trans,ekin_rot,ekin_tot,jumpH"/>
  <raster
      metric="ekin_tot"
      resolution="2.0"
      statistic="percentile"
      p="95"/>
  <raster
      metric="jumpH"
      resolution="2.0"
      statistic="percentile"
      p="95"/>
  <raster
      metric="count"
      resolution="2.0"/>
  <shape metric="deposition"/>
</output>

<debug
    timing_report="true"/>
```

## Outputs

### Text output

If no `<output>...</output>` tag is given, the simulator writes text trajectory files and a `manifest.csv`.

### Binary output

If `<output><trajectory format="binary" .../></output>` is given, the simulator writes:

- `manifest.csv`
- one binary shard file per worker thread, for example `trajectories_001.rfbin`

The binary layout is schema-driven and depends on the trajectory `metric` list in the XML.

### Raster output

If one or more `<raster .../>` tags are given inside `<output>`, the simulator writes one GeoTIFF per raster definition.

Examples:

- `q95_ekin_tot.tif`
- `q95_jump_h.tif`
- `count_rocks.tif`

### Deposition points

If `<shape metric="deposition" .../>` is given inside `<output>`, the simulator writes one point shapefile, for example:

- `deposition_points.shp`

Each point stores the final accepted trajectory coordinate as 3D point geometry and writes additional attributes into the DBF.

### Timing report

If `<debug timing_report="true"/>` is enabled, the simulator also writes:

- `timing_report.txt`

This contains batch-level timings and peak working-set memory.

## References

Conceptual references relevant to the implemented methods are listed in [`REFERENCES.md`](REFERENCES.md).

## Third-Party Code and Attribution

Bundled third-party source code currently found in this repository:

- Quickhull by Anatoliy V. Tomilov

See [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the exact required notice text for bundled third-party code.

External libraries used by the build but not bundled in this repository include:

- Eigen
- GDAL

## License

Add your project license here before public release.
