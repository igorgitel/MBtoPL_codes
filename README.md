# MBtoPL_codes

This repository provides the full simulation suite used in the MBtoPL study:  
**"When does Maxwell-Boltzmann transforms to a Power-Law?"**

Each folder contains a distinct physical scenario that explores how relativistic collisions, injections, and cascades lead to the formation of nonthermal energy distributions from initially thermal (or monoenergetic) populations.


---

## 📁 Folders and Simulation Scenarios

### `single_mass_close_box/`
A closed system of identical particles undergoing elastic collisions. Serves as the baseline thermalization test toward a Maxwell–Jüttner distribution.

### `two_type_mass_closed_box/`
Simulates a two-mass system with no injection. Demonstrates energy exchange dynamics between light and heavy particles through relativistic collisions.

### `Asceleration_open_box/`
Simulates stochastic acceleration of particles in an open system. A single bunch of low-mass particles is injected into a static background of very high-mass particles. The background particles do not evolve — they serve as an energy reservoir with fixed velocities and directions. The code tracks how the low-mass particles gain energy through repeated elastic collisions with this background. This scenario is designed to model non-relativistic energy exchange and tests the distribution function of this low-mass particles.

### `Relativistic_Asceleration_open_box/`
Simulates stochastic acceleration of particles in an open system, similar to the `Asceleration_open_box/` setup but using fully relativistic kinematics. A single bunch of low-mass particles is injected into a background of static, high-mass particles that do not evolve. These background particles act as a fixed-energy reservoir, and only the low-mass population is tracked. The code models how energy is transferred to the low-mass particles through repeated relativistic elastic collisions, allowing analysis of how relativistic energy drift shapes their distribution function.


### `Ultra_Relativistic_Asceleration_on_infinite_mass_open_box/`
Simulates stochastic acceleration in the **ultra-relativistic limit**. A bunch of finite-mass particles interacts with a background of **infinite-mass scatterers**, all moving at a fixed relativistic velocity. The scatterers do not change direction or energy, and no recoil is applied — this isolates the effect of energy transfer on the evolving particle population. The code follows how finite-mass particles gain energy through repeated **relativistic collisions** with this single-velocity magnitude, rigid background.


### `Cascade_relativistic_dynamics_open_box/`

Implements a **single monoenergetic injection** of particles into a thermal background, designed to model energy redistribution via relativistic elastic collisions. The injected particles interact with the background, gradually transferring energy through collisions. As the system evolves, more background particles become energetically active, leading to a broadening of the energy distribution — a process analogous to a **cascade**.


### `Energy_injection_to_thermal_Kappa_open_box/`

This code simulates a **open system** of colliding particles with a **fixed total number**, where energy is periodically injected into a random subset of the particles. At regular time intervals, a group of particles is selected and their energies are reassigned from a high-temperature thermal distribution.

The entire system evolves through relativistic elastic collisions. The code tracks how the energy distribution evolves over time and examines whether a steady state with a nonthermal tail—such as a kappa distribution—might emerge as a result of the continuous energy input, although this is not necessarily the case.

---

## 🧠 Physics Goals

- Understand under what conditions relativistic systems transition from Maxwellian to power-law distributions
- Compare equilibrium (closed-box) vs. driven (open-box or injected) scenarios
- Explore the role of mass ratios, injection energy, and cross-section dependence

---

## 📌 Dependencies

All folders use shared utilities located in `Basic functions/`, including:
- `collision.m`
- `collision_freq.m`
- `histlog.m`
- `randdir_matrix.m`

---


## 🧪 How to Run

Each folder contains a self-contained script or a named main function.
Open Matlab and run witin the folder. Make sure all the files witin the folder are present.
