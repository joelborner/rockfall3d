# rockfall3d

`rockfall3d` is a C++ rigid-body rockfall simulator for batch trajectory simulations on raster DEMs with roughness.

The public-facing workflow is XML-driven:

- choose a DEM raster
- choose one or more rock point clouds
- choose release points
- optionally define roughness
- choose solver (JOR-prox, SOR-prox)
- choose number of parallel threads
- write trajectories in text or binary form

## Features

- DEM-based terrain interaction with bilinear interpolation
- rigid-body dynamics with quaternion orientation
- JOR/SOR-prox contact solve for terrain interaction
- roughness model
- deterministic orientation and roughness seeds
- shared-memory, thread-based parallel execution
- text or binary trajectory output
- optional terrain export and debugging options

## Repository Layout

- [`rockfall3d/`](rockfall3d): source code
- [`sample_data/`](sample_data): example XML files and sample inputs
- [`rockfall3d.sln`](rockfall3d.sln): Visual Studio solution
- [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md): bundled third-party attribution
- [`REFERENCES.md`](REFERENCES.md): conceptual references for important methods used in the code

## Dependencies

### Bundled in the repository

- Quickhull by Anatoliy V. Tomilov, bundled in [`rockfall3d/external/quickhull/include/quickhull.hpp`](rockfall3d/external/quickhull/include/quickhull.hpp)

### External dependencies

- Visual Studio 2022 with MSVC toolset (`v143`)
- Eigen 3.4.x
- GDAL
- vcpkg is recommended for installing GDAL on Windows

Notes:

- The checked-in Visual Studio project currently builds with `WITH_GDAL` enabled.
- Quickhull is bundled, so it does not need to be installed separately.
- The project file currently references Eigen through a local include path; you may need to update that path on your machine.

## Installation on Windows

### 1. Install Visual Studio

Install Visual Studio 2022 with the Desktop development with C++ workload.

### 2. Install vcpkg and GDAL

Install [vcpkg](https://github.com/microsoft/vcpkg), then install GDAL for x64 Windows:

```powershell
vcpkg install gdal:x64-windows
vcpkg integrate install
```

### 3. Install Eigen

Download Eigen 3.4.x and make sure the include directory is available to the project.

The current project file references Eigen via `AdditionalIncludeDirectories` in [`rockfall3d/rockfall3d.vcxproj`](rockfall3d/rockfall3d.vcxproj), so you may need to update that path to match your local installation.

### 4. Build

From Visual Studio:

1. Open [`rockfall3d.sln`](rockfall3d.sln)
2. Select `Release | x64`
3. Build the solution

The executable is written to:

```text
x64\Release\rockfall3d.exe
```

## Running

The executable expects exactly one argument: the XML config file.

```powershell
.\rockfall3d.exe path\to\config.xml
```

## Input Files Referenced by XML

### Rock list file

The file referenced by `<rocks list="..."/>` contains one rock point-cloud path (`.pts`) per line.

- relative paths are resolved relative to the list file location

### Start points file

The file referenced by `<starts list="..."/>` contains:

```text
X Y dz
```

where:

- `X`, `Y` are plan coordinates
- `dz` is the desired initial clearance between the terrain and the lowest rock vertex after the initial orientation is applied

## XML Configuration Reference

The parser expects self-closing tags.

Required tags:

- `<terrain/>`
- `<rocks/>`
- `<starts/>`
- `<simulation/>`
- `<stop/>`
- `<solver/>`

Optional tags:

- `<roughness/>`
- `<output/>`
- `<debug/>`

The top-level `<rockfall3d/>` marker is accepted in sample files, but the parser does not rely on it.

### `<terrain/>`

Required attributes:

| Attribute | Meaning |
|---|---|
| `dem` | Path to the DEM raster (`.tif`, `.tiff`, or `.asc`) |

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
| `mu` | Coulomb friction coefficient |
| `dt` | Time step |
| `output_dir` | Output directory |

Optional attributes:

| Attribute | Default if omitted | Meaning |
|---|---|---|
| `dump_step` | falls back to `dt` | Output sampling interval |
| `parallel_threads` | `1` | Number of parallel worker threads requested |
| `random_seed` | random device | Fallback seed for both orientation and roughness if dedicated seeds are omitted |
| `orientation_seed` | `random_seed` if provided, otherwise random device | Seed for random orientations |
| `roughness_seed` | `random_seed` if provided, otherwise random device | Seed for roughness realization |

Notes:

- If `parallel_threads="0"` is given, it is clamped to `1`.
- The code clamps the requested thread count to the number of available logical processors.

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

### `<output/>`

This tag is optional.

If omitted:

- output falls back to legacy text trajectory files

If present:

Required attributes when tag is present:

| Attribute | Meaning |
|---|---|
| `format` | Output format: `binary-sharded` or `legacy-text` |
| `variables` | Comma-separated output variables |

Valid variable names / aliases:

| Variable | Meaning |
|---|---|
| `t` | time [s] |
| `x` | X coordinate [m] |
| `y` | Y coordinate [m]|
| `z` | Z coordinate [m]|
| `quat` | Orientation (Quaternion) |
| `vx,vy,vz` | Translational velocities [m/s] |
| `vtot` | Resultant translational velocity [m/s]  |
| `wx,wy,wz` | Rotational velocities [rad/s]|
| `wtot` | Resultant rotational velocity [rad/s] |
| `ekin_trans` | Kinetic translational energy [J] |
| `ekin_rot` | Kinetic rotational energy [J] |
| `ekin_tot` | Kinetic energy [J] |
| `jumpH` | Jump height [m] |

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

<terrain dem="data/dem_2m.tif"/>

<rocks
    list="data/rocks.txt"
    density="2500"
    random_orientations="100"/>

<starts list="data/start_points.txt"/>

<simulation
    mu="0.5"
    dt="0.002"
    dump_step="0.02"
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

<solver
    type="jor-prox"
    alpha="0.5"
    atol="1e-8"
    rtol="1e-6"
    maxiter="10000"/>

<roughness
    raster="data/roughness_grid_size.tif"
    step_func="smoothstep"/>

<output
    format="binary-sharded"
    variables="time,x,y,z,qw,qx,qy,qz,vx,vy,vz,vtot,wx,wy,wz,wtot,ekin_trans,ekin_rot,ekin_tot,jumpH"/>

<debug
    timing_report="true"/>
```

## Outputs

### Legacy text output

If no `<output/>` tag is given, the simulator writes legacy text trajectory files and a `manifest.csv`.

### Binary-sharded output

If `<output format="binary-sharded" .../>` is given, the simulator writes:

- `manifest.csv`
- one binary shard file per worker thread, for example `trajectories_001.rfbin`

The binary layout is schema-driven and depends on the `variables` list in the XML.

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
