# MBtoPL_codes

This repository contains MATLAB simulation tools supporting the MBtoPL project—**"From Maxwell–Boltzmann to Power Law: Energy Cascades in Relativistic Particle Distributions"**. The code is designed to explore how monoenergetic or thermal (“Maxwell–Boltzmann”) particle injections evolve into non-thermal power-law spectra through stochastic relativistic collisions.

---

## 📚 Connection to MBtoPL Paper

The simulations here directly implement the mathematical model described in the MBtoPL paper:

- **`func_One_type_prep.m`**, **`func_asc_single_injection_Mat_TS.m`**, and **`SPI_RED_cascade.m`** simulate frameworks of monoenergetic and continuous injection relaxation in relativistic plasmas.
- They reproduce the core findings of the paper:
  - injection of delta-function or thermal distributions into a background
  - subsequent evolution into Maxwell–Jüttner or extended power-law tails
  - effects of cross-section energy dependence and injection strength
- The code corroborates analytic and semi-analytic MBtoPL predictions and produces plots comparable to figures in the paper.
