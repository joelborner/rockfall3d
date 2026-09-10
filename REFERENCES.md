# References

This file lists conceptual references relevant to methods used in `rockfall3d`.

These references are included as background for important implemented ideas. They are not claims that code was copied from these publications unless explicitly stated in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Gradient Noise / Perlin Roughness

- Perlin, K. (1985). *An Image Synthesizer*. ACM SIGGRAPH Computer Graphics, 19(3), 287-296.

Relevant in this project for:

- roughness built from globally defined Perlin fields

## Quickhull Concept

- Barber, C. B., Dobkin, D. P., & Huhdanpaa, H. (1996). *The Quickhull Algorithm for Convex Hulls*. ACM Transactions on Mathematical Software, 22(4), 469-483.
- Mehlhorn, K., Naeher, S., Schilz, T., Schirra, S., Seel, M., Seidel, R., & Uhrig, C. (1996). *Checking geometric programs or verification of geometric structures*. In Proceedings of the 12th Annual ACM Symposium on Computational Geometry, 159-165.

Relevant in this project for:

- convex hull construction of rock point clouds prior to mass/inertia computation

Note:

- the repository bundles a specific Quickhull implementation by Anatoliy V. Tomilov; the exact notice text for that bundled source is listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
- the Tomilov Quickhull repository also points to the Mehlhorn et al. reference above in the context of geometric verification

## Polyhedral Mass Properties

- Mirtich, B. (1996). *Fast and Accurate Computation of Polyhedral Mass Properties*. Journal of Graphics Tools, 1(2), 31-50.

Relevant in this project for:

- mass, center-of-mass, and inertia computation from the rock convex hull using tetrahedral accumulation ideas

## Nonsmooth Rigid-Body Dynamics and Hard Contact

- Acary, V., & Bourrier, F. (2021). *Coulomb friction with rolling resistance as a cone complementarity problem*. European Journal of Mechanics - A/Solids, 85, 104046. https://doi.org/10.1016/j.euromechsol.2020.104046
- Leine, R. I., Schweizer, A., Christen, M., Glover, J., Bartelt, P., & Gerber, W. (2014). *Simulation of rockfall trajectories with consideration of rock shape*. Multibody System Dynamics, 32, 241-271. https://doi.org/10.1007/s11044-013-9393-4
- Leine, R. I., & Nijmeijer, H. (2004). *Dynamics and Bifurcations of Non-Smooth Mechanical Systems*. Springer.
- Leine, R. I., & van de Wouw, N. (2008). *Stability and Convergence of Mechanical Systems with Unilateral Constraints*. Springer.
- Leine, R. I. *An Introduction to Nonsmooth Dynamics*, to be published by Springer.

Relevant in this project for:

- rigid-body time stepping with unilateral contact and frictional impact/contact laws
- the slippage-dependent Coulomb coefficient and its internal-state evolution
- the per-contact rolling-friction cone, maximum-dissipation law, and closed-form cone projection
- the general nonsmooth dynamics framework behind the contact formulation
