% Define simulation metadata to locate the file
N = 1e5;             % Number of Type 1 particles
icb.beta = 0.0995;   % Velocity (v/c) of background scatterers
ic.Energy0 = 100;    % Initial energy of particles (e.g., gamma * m = 100)
time_step = 1.0;     % Time between histogram snapshots
max_time = 1e4;      % Total simulation time (in code units)

data = Single_injection_evolution_ulta_rel(N,ic,icb,max_time,time_step);
