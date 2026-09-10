# Inclined-plane cube test

This minimal case tests whether a non-spherical rigid body can lose contact
with a geometrically smooth plane using constant Coulomb friction only.

- Terrain: LV95 (`EPSG:2056`) GeoTIFF containing an exact 45-degree plane,
  descending in positive X
- Rock: homogeneous 1 m cube at 2500 kg/m3
- Initial state: zero velocity, deterministic random orientation, 1 mm clearance
- Contact: constant Coulomb coefficient `mu = 0.8`
- Restitution: `en = et = 0`
- Excluded: roughness, rolling resistance, and slippage-dependent friction

Run from the repository root:

```powershell
.\build\windows-release\Release\rockfall3d.exe .\sample_data\inclined_plane_cube\inclined_plane_cube.xml
```

The text trajectory is written to `output`. Positive `jumph` values after the
first terrain contact indicate separation. Because the cube starts in a
non-face-aligned orientation, this case tests geometry-induced tumbling rather
than the special symmetric case of a cube resting exactly on one face.
