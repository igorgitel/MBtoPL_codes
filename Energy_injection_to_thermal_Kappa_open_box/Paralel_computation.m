clear
clc
addpath('Basic functions/')
delete(gcp('nocreate')) % close existing pool
parpool(9);             % create new pool with 5 workers

%% Simulation parameters
max_time = 1e2;                     % total simulation time
cs_p = 0;                           % power-law index for cross section: sigma ~ E^(-cs_p)
for time_between_injections = 10.^(-3:1)          % dt  % time step between energy injections
num_of_inj_particles = 10;           % n % Number of particles to receive energy per injection
energy_inc = 1;                     % e % energy increment during each injection
num_runs = 18;                     % number of simulation repetitions
is_save=false;                      % Whether to save each individual run
save_folder='';                     % Folder to save results (unused if is_save = false)
parallel_computation = true; 

 
%% Particle initial conditions
ic.m = 1;                             % particle mass
ic.N = 1e7;                           % number of particles
ic.gb = 0.05;                         % beta = v/c
ic.g = sqrt(ic.gb^2 + 1);             % Lorentz factor
ic.e = ic.m * ic.g;                   % total energy
ic.ek = ic.e - ic.m;                  % kinetic energy
ic.edges = 10.^(-16:0.1:16);          % energy bin edges (logarithmic)

Heat_rate = energy_inc * num_of_inj_particles / time_between_injections;
% disp(['Heating rate = ' num2str(Heat_rate) ' energy units per time unit'])
heat_rate = Heat_rate / ic.N;  % heating rate per particle
  

%% Calculate Maxwell-Jüttner temperature T
syms x
T = double(vpasolve(3*x + ic.m * besselk(1,ic.m/x) / besselk(2,ic.m/x) == ic.e, x));
clear x

%% Prepare cell array to store results from all runs
all_data = cell(1, num_runs);

%% Run all simulations in parallel
parfor k = 1:num_runs
    tic
    rng(k);  % Set independent random seed for each run
    disp(['run ' num2str(k)])
    ic_k = ic;  % make a local copy for each worker
    all_data{k} = Energy_injection_no_cooling(ic, ...
                                           cs_p, ...
                                           energy_inc, ...
                                           num_of_inj_particles,...
                                           time_between_injections, ...
                                           max_time, ...
                                           save_folder, ...
                                           is_save, ...
                                           parallel_computation);
    toc
end

%% Assume all runs have the same number of time steps
num_times = length(all_data{1});
f_sim_sum = cell(1, num_times);         % sum of f_sim over all runs
bins_ref = all_data{1}(1).bins;         % reference energy binning

% Initialize accumulation array
for i = 1:num_times
    f_sim_sum{i} = zeros(size(all_data{1}(i).f_sim));
end

% Sum f_sim across all runs for each time step
for k = 1:num_runs
    for i = 1:num_times
        f_sim_sum{i} = f_sim_sum{i} + all_data{k}(i).f_sim;
    end
end

% Average f_sim for each time step
f_sim_avg = cell(1, num_times);
data = all_data{1};  % use structure from the first run as template

for i = 1:num_times
    data(i).f_sim = f_sim_sum{i} / num_runs;
end

%% Save averaged result with all relevant parameters
save_folder = ['mean_distribution_from_' num2str(num_runs) '_runs_3'];
if ~exist(save_folder, 'dir')
    mkdir(save_folder);  % create folder if it doesn't exist
end
% load(fullfile(save_folder, 'mean_f_sim_data.mat'))
str =  ['N_' num2str(ic.N,'%1.0e') ...
       '_e_' num2str(energy_inc,'%1.0e') ...
       '_dt_' num2str(time_between_injections,'%1.0e') ...
       '_n_' num2str(num_of_inj_particles,'%1.0e') ...
       '_H_' num2str(Heat_rate,'%1.1e')...
       '_max_time_' num2str(max_time,'%1.0e') ...
       '_cs_p_' num2str(cs_p)];

save(fullfile(save_folder, [str '.mat']), ...
    'data', 'ic', 'cs_p', 'energy_inc', ...
    'time_between_injections', 'max_time', ...
    'bins_ref', 'T','num_of_inj_particles', ...
    'Heat_rate');
end
%% Animation script 
% anima_bach_render