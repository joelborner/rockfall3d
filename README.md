# rockfall3d

`rockfall3d` is a C++ rigid-body rockfall simulator for batch trajectory simulations on raster DEMs with roughness.

## Contact
Joël Borner

WSL Institute for Snow and Avalanche Research SLF

joel.borner@slf.ch

## Features

- DEM-based terrain interaction with bilinear interpolation
- rigid-body dynamics with quaternion orientation
- nonsmooth dynamics with hard contact laws and Coulomb friction. Moreau timestepping scheme and JOR/SOR-prox contact solve for terrain interaction
- roughness model based on Perlin noise
- deterministic random orientation and roughness, reproducible with seeds
- shared-memory, thread-based parallel execution
- trajectory and/or raster output
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
X Y Z dz
```

where `dz` is optional:

- if 3 columns are defined, `X`, `Y`, `Z` are the absolute start coordinates of the centre of mass of the rock.
- if 4 columns are defined,`dz` is the desired initial clearance between the terrain and the lowest rock vertex after the initial orientation is applied. In this case, `Z` is ignored. 

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
- `<output>...</output>`
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

Each point stores the final accepted trajectory coordinate as 3D point geometry and writes the following attributes into the DBF:

- `run_idx`
- `rock_idx`
- `rock_path`
- `start_idx`
- `orient_idx`
- `start_x`
- `start_y`
- `start_z`
- `status`
- `stop_crit`

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
