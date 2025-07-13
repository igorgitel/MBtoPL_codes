func_One_type_prep.m

This script performs a relativistic Monte Carlo simulation of colliding identical particles in a closed box with a single mass species. It models the time evolution of particle energies under two-particle elastic collisions, starting from a uniform initial energy distribution.

🔧 Inputs:
	•	m — particle rest mass
	•	gb — initial dimensionless bulk velocity (used to compute Lorentz factor)
	•	N — number of particles in the simulation
	•	step_inc — time step interval for data sampling

📤 Outputs:
	•	Enb — final total energy (including rest mass) of each particle at end of simulation
	•	data — struct array containing time-series snapshots of:
	•	mean kinetic energy
	•	energy and momentum histograms
	•	theoretical distributions (relativistic Maxwell–Jüttner)
	•	collision count
	•	simulation time
	•	icb — struct with simulation constants (mass, initial energy, binning info)

⚙️ Core Features:
	•	Initializes particles with the same total energy derived from a Lorentz-boosted rest energy
	•	Computes theoretical temperature T by solving a transcendental equation involving Bessel functions
	•	Uses random isotropic directions for momentum sampling (randdir_matrix)
	•	Evolves the system via stochastic binary collisions using the RelCol function
	•	Uses an adaptive time loop with a probabilistic collision acceptance step
	•	Samples kinetic energy and momentum (μ) distributions at regular intervals
	•	Compares simulation histograms with the Maxwell–Jüttner theoretical distribution

📦 Output:
	•	Simulation results are saved in a .mat file in the saves/ directory, with a filename summarizing parameters

⸻

