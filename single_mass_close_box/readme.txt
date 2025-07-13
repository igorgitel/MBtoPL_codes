func_One_type_prep.m

This script performs a relativistic Monte Carlo simulation of colliding identical particles in a closed box with a single mass species. It models the time evolution of particle energies under two-particle elastic collisions, starting from a uniform initial energy distribution.

🔧 Inputs:
	•	m — particle rest mass
	•	gb — initial bulk velocity in units of c (used to compute initial Lorentz boost)
	•	N — number of particles in the simulation
	•	step_inc — time step between data samples

⸻

📤 Outputs:
	•	Enb — vector of total energies (rest + kinetic) for all particles at the final simulation time
	•	data — structure array of sampled time snapshots with:
	•	kin_energy: mean kinetic energy per particle
	•	bins_Ek, f_sim_Ek: energy histogram and distribution
	•	Teo_Ek: Maxwell–Jüttner distribution for comparison
	•	bins_mu, f_sim_mu: distribution of particle momentum (μ)
	•	Teo_mu: theoretical prediction for μ
	•	coll_num: cumulative number of accepted collisions
	•	time: time of snapshot
	•	icb — struct with simulation parameters (mass, energy, bin edges)

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

% Example Usage in MATLAB

% Set rest mass of particles
m = 1;

% Set initial scaled momentum (dimensionless, v/c)
gb = 0.1;

% Convert to Lorentz boost factor (gamma)
b = sqrt(gb^2 / (1 + gb^2));
g = 1 / sqrt(1 - b^2);

% Calculate initial energy per particle (total and kinetic)
e = m * g;
ke = e - m;  % kinetic energy

% Number of particles
N = 1e6;

% Time step interval for data collection
step_inc = 1e-3;

% Run the simulation
[Enb, data, icb] = func_One_type_prep(m, gb, N, step_inc);
