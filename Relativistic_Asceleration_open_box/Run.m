clear 
clc 
close all
%% 🧊 STEP 1: Prepare background (high-mass thermal bath)

icb.m  = 1e8;          % 🪨 Background particle mass
icb.N  = 1e6;           % 🔢 Number of background particles
icb.gb = 1e-1;          % 🌀 Initial reduced momentum (gamma * beta)

step_inc = 1e-1;        % 🧭 Time resolution for saving diagnostics (background prep)
max_time = 6;           % ⏱️ Duration of background evolution

% ⚗️ Generate and save background energy distribution (thermal MJ)
% This initializes the thermal bath used later in the injection simulation

str = ['One_type_mass_' num2str(icb.m, '%1.0e') ...
       '_initial_gb_' num2str(icb.gb, '%1.0e') ...
       '_particle_number_' num2str(icb.N, '%1.0e')];

file_path = fullfile('saves_one_type', [str '.mat']);

if exist(file_path, 'file')
    disp('File exists.');
else 
    warning('File not found!');
    [~, ~, ~] = func_One_type_prep(icb.m, icb.gb, icb.N, step_inc, max_time);
end
%% 🟡 STEP 2: Define initial condition for the injected (low-mass) population

% The injected particles are initialized with a delta-function distribution in energy
ic.gb = 1e-1;                                % 🚀 Reduced momentum of injected particles
ic.b  = sqrt(ic.gb^2 / (1 + ic.gb^2));       % 🔄 Velocity derived from gb
ic.g  = 1 / sqrt(1 - ic.b^2);                % 🔁 Lorentz gamma factor

ic.m  = 1;                                   % 🌱 Mass of injected particles
ic.e  = ic.g * ic.m;                         % ⚡ Total energy per particle
ic.ke = ic.e - ic.m;                         % 🔥 Kinetic energy

%% ⚙️ Simulation Parameters

f_param     = 1e1;     % 🔧 Scaling factor for collision frequency
num_par_inj = 1e6;     % 💡 Number of injected low-mass particles
max_time    = 1e4;     % ⏱️ Duration of the injection simulation
step_inc    = 1.0;     % 🧭 Time between saved snapshots

%% 🔁 STEP 3: Run the simulation — single monoenergetic injection

folder = 'saves_v1';   % 💾 Output folder for saving results

[data, En] = func_asc_single_injection_Mat_TS(...
    ic, ...            % Injected particle parameters
    icb, ...           % Background particle parameters
    f_param, ...       % Collision frequency scaling
    num_par_inj, ...   % Number of injected particles
    max_time, ...      % Simulation time
    step_inc, ...      % Time between saved snapshots
    folder ...         % Save path
);
