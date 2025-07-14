% Add path to folder containing helper functions
addpath('Basic functions/')

%% 🔧 Background (Type 2) Particle Setup

% Rest mass of background particles
icb.m = 1;

% Initial reduced momentum (gamma * beta) of background particles
icb.gb = 1e-4;

% Number of background particles
icb.N = 1e5;

% Compute Lorentz gamma from reduced momentum: γ = sqrt(1 + (p/m)^2)
icb.g = sqrt(icb.gb^2 + 1);

% Total energy per background particle: E = γ * m
icb.e = icb.m * icb.g;

% Kinetic energy per particle: KE = E - m
icb.ke = icb.e - icb.m;

% Folder where the background particle distribution will be saved
icb.folder = 'One_type';

% Generate and save background energy distribution using func_One_type_prep
[Enb, ~, ~] = func_One_type_prep(icb.m, icb.gb, icb.N, icb.folder);

%% 🔷 Injected Particle Setup (Type 1)

% Initial reduced momentum (gamma * beta) of injected particles
ic.gb = 1e7;          

% Injected particle mass (equal to background for this setup)
ic.m = icb.m;

%% ⏱ Simulation Parameters

max_time     = 1e4;    % Total simulation time (arbitrary units)
maxnumofpar  = 1e5;    % Maximum total number of injected particles allowed
numofpar     = 1e4;    % Number of particles injected at each injection step

% Folder to save simulation results
folder = 'RD_saves_05';

%% 🚀 Run the Simulation (Relativistic Energy Drift Cascade)

tic
[data, str, ic, icb] = SPI_RED_cascade(ic, icb, ...
    maxnumofpar, numofpar, ...
    max_time, folder);
toc

%% 📊 Optional: Animate spectral evolution
% evol_Mat_anim(data, str, folder, ic, icb)
