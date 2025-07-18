%----------------------------------------------------------
% Simulation of two-species relativistic collisions
% using par_col.m to track thermalization and energy sharing
%----------------------------------------------------------

gb = 1e-3;             % Reduced momentum (gamma * beta) for both species

fparam = 1e1;          % Scaling factor for collision frequency (dimensionless)

% --- Define properties of the first particle species (light particles) ---
ic1.m = 1;             % Mass of species 1
ic1.gb = gb;           % Initial reduced momentum of species 1

% --- Define properties of the second particle species (heavy particles) ---
ic2.m = 1836;          % Mass of species 2 (e.g., proton/electron mass ratio)
ic2.gb = gb;           % Same reduced momentum (initial velocity) as species 1

N = 1e6;               % Number of particles per species (total = 2N)
max_time = 1e4;        % Total simulation time (in arbitrary units)

% --- Run simulation ---
% This function models a two-species system with relativistic elastic collisions.
% Each species starts monoenergetic; the simulation tracks their energy evolution.
% Outputs:
% - data: struct with time history of energy distributions and collision counts
% - ic1, ic2: updated structures for each species, including computed energies
[data, ic1, ic2] = par_col(ic1, ic2, N, max_time, fparam);
