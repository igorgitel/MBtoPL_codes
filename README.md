# MBtoPL_codes

This repository contains MATLAB simulation tools supporting the MBtoPL project—**"From Maxwell–Boltzmann to Power Law: Energy Cascades in Relativistic Particle Distributions"**. The code is designed to explore how monoenergetic or thermal (“Maxwell–Boltzmann”) particle injections evolve into non-thermal power-law spectra through stochastic relativistic collisions.

---

## 📚 Connection to MBtoPL Paper

This repository is organized by physical scenario, with each folder or script corresponding to a distinct kinetic setup:

- **`single_mass_close_box/`** simulates an idealized **closed-box gas** of identical particles undergoing random relativistic collisions. It models thermalization within a fixed particle pool and serves as the baseline case for equilibrium relaxation toward the Maxwell–Jüttner distribution.

- **`TwotypeColision.m`** models **stochastic acceleration** of particles via repeated elastic collisions with a moving background. This includes both **non-relativistic and relativistic target velocities**, and tracks how initially cold particles gain energy through interactions with a prescribed population of scatterers.

- **`SPI_RED_cascade.m`** explores **cascade-like energy redistribution**: a monoenergetic population is injected once into a thermal background, and energy is spread via repeated collisions. The system begins in a **mildly relativistic state** and naturally evolves toward a broader energy distribution, mimicking realistic plasma heating scenarios.

- **`One_type_prep_with_energy_change.m`** implements a scenario where energy is **continuously injected into a fixed number of particles** at specified intervals. This enables controlled, local energy driving and may reproduce **nonthermal steady-state distributions**, including **kappa-like tails**, depending on the injection frequency and strength.

Each module isolates a different aspect of kinetic behavior—from idealized equilibration to driven turbulence and collisional energization—allowing targeted investigation of the conditions under which a Maxwellian distribution transitions to a power law.
