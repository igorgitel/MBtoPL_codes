clear
clc
close
addpath('Basic functions/')           % Add folder with custom functions to the MATLAB path

%% Simulation parameters
max_time = 2e2;                       % Total simulation time
cs_p = 0;                             % Power-law index for the cross section: sigma ∝ E^(-cs_p)

time_between_injections = 1           % Time interval between energy injections
energy_inc = 1;                       % Energy injected per event
num_of_inj_particles = 1;             % Number of particles to receive energy per injection

Heat_rate = energy_inc * num_of_inj_particles / time_between_injections;
% disp(['Heating rate = ' num2str(Heat_rate) ' energy units per time unit'])

%% Particle initial conditions
ic.m = 1;                             % Particle mass
ic.N = 1e6;                           % Total number of particles
ic.gb = 0.05;                         % Initial velocity in units of c (v/c)
ic.g = sqrt(ic.gb^2 + 1);             % Lorentz factor
ic.e = ic.m * ic.g;                   % Total energy
ic.ek = ic.e - ic.m;                  % Initial kinetic energy
ic.edges = 10.^(-16:0.1:16);          % Energy bin edges for histogram (logarithmic scale)

folder = 'Background_population';
str = ['MJ_Mass_' num2str(ic.m,'%1.0e') ...
       '_gb_' num2str(ic.gb,'%1.0e') ...
       '_partcle_number_' num2str(ic.N,'%1.0e')];
file_path = fullfile(folder, [str '.mat']);

% Убедимся, что папка существует
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Проверка наличия файла
if ~(exist(file_path, 'file')==2)
    % Файл не существует — создаём его
    Enb = func_One_type_prep(ic.m, ic.gb, ic.N);
else
    % Файл уже существует — уведомляем
    [~, fname, ext] = fileparts(file_path);
    disp(['File "' fname ext '" already exists in folder "' folder '".'])
end

%% Output configuration
save_folder = 'saves_energy_inj_no_cooling_v2';  % Folder where results will be saved
is_save = true;                            % Set to true to save the output to a .mat file
parallel_computation =false;
%% Run the simulation (single run)
data = Energy_injection_no_cooling(ic, ...
                                  cs_p, ...
                                  energy_inc, ...
                                  num_of_inj_particles, ...
                                  time_between_injections, ...
                                  max_time, ...
                                  save_folder, ...
                                  is_save, ...
                                  parallel_computation);
%% 
% str =  ['N_' num2str(ic.N,'%1.0e') ...
%        '_e_' num2str(energy_inc,'%1.0e') ...
%        '_dt_' num2str(time_between_injections,'%1.0e') ...
%        '_n_' num2str(num_of_inj_particles,'%1.0e') ...
%        '_H_' num2str(Heat_rate,'%1.1e')...
%        '_max_time_' num2str(max_time,'%1.0e') ...
%        '_cs_p_' num2str(cs_p)];

% anima_s(data,ic,Heat_rate,str,save_folder)


anima_bach_render
