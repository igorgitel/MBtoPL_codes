
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
Models stochastic acceleration of particles due to structured collisions in an open system. Designed to test gradual energy gain via random scattering.

### `Relativistic_Asceleration_open_box/`
Tests scenarios where the ascelerating population is relativistic. Demonstrates how relativistic energy contributes to nonthermal energy growth.

### `Ultra_Relativistic_Asceleration_on_infinite_mass_open_box/`
Simulates acceleration of finite-mass particles interacting with a background of infinite-mass scatterers.

### `Cascade_relativistic_dynamics_open_box/`
Implements a monoenergetic injection into a thermal background. The resulting cascade models how injected energy spreads.

### `Injection_of_high_energy_termal_kappa_close_box/`
Simulates repeated energy injections into a closed system. Aimed at reproducing steady-state kappa distributions through continuous driving and collisions.

---

## 🧠 Physics Goals

- Understand under what conditions relativistic systems transition from Maxwellian to power-law distributions
- Compare equilibrium (closed-box) vs. driven (open-box or injected) scenarios
- Explore the role of mass ratios, injection energy, and cross-section dependence

---

## 🧪 How to Run

Each folder contains a self-contained script or a named main function.
Open Matlab and run witin the folder. Make sure all the files witin the folder are present.
