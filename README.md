
# MBtoPL_codes

This repository provides the full simulation suite used in the MBtoPL study:  
**"From Maxwell–Boltzmann to Power-Law: Energy Cascades in Relativistic Particle Distributions."**

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
Tests scenarios where both colliding populations are relativistic. Demonstrates how relativistic energy drift contributes to nonthermal energy growth.

### `Ultra_Relativistic_Asceleration_on_infinite_mass_open_box/`
Simulates acceleration of finite-mass particles interacting with a background of infinite-mass scatterers. Used to isolate pure energy gain without recoil.

### `Cascade_relativistic_dynamics_open_box/`
Implements a monoenergetic injection into a thermal background. The resulting cascade models how injected energy spreads and leads to power-law tails.

### `Injection_of_high_energy_termal_kappa_close_box/`
Simulates repeated energy injections into a closed system. Aimed at reproducing steady-state kappa distributions through continuous driving and collisions.

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

Each folder contains a self-contained script (usually `run.m`) or a named main function. You can:
1. `cd` into the folder
2. Open MATLAB
3. Run the corresponding script
4. Use built-in plotting commands or save data for comparison to the MBtoPL figures

---

For questions, extensions, or collaborations, please reach out or open an issue.
