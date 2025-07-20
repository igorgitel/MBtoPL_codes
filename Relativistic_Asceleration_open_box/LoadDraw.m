clear
close all
clc

% background (high-mass thermal bath)
icb.m  = 1e10;          % 🪨 Background particle mass
icb.N  = 1e6;           % 🔢 Number of background particles
icb.gb = 1e-3;          % 🌀 Initial reduced momentum (gamma * beta)

% The injected particles are initialized with a delta-function distribution in energy
ic.gb = 1e-3;                                % 🚀 Reduced momentum of injected particles
ic.b  = sqrt(ic.gb^2 / (1 + ic.gb^2));       % 🔄 Velocity derived from gb
ic.g  = 1 / sqrt(1 - ic.b^2);                % 🔁 Lorentz gamma factor

ic.m  = 1;                                   % 🌱 Mass of injected particles
ic.e  = ic.g * ic.m;                         % ⚡ Total energy per particle
ic.ke = ic.e - ic.m;                         % 🔥 Kinetic energy

% ⚙️ Simulation Parameters

f_param     = 1e1;     % 🔧 Scaling factor for collision frequency
num_par_inj = 1e6;     % 💡 Number of injected low-mass particles
max_time    = 1e4;     % ⏱️ Duration of the injection simulation
step_inc    = 1.0;     % 🧭 Time between saved snapshots

data_save_folder = 'saves_v1';   % 💾 Output folder for saving results

str=['mass_Background_' num2str(icb.m,'%1.0e')...
    '_Single_injection_f_param_' num2str(f_param,'%1.0e')...
    '_num_par_inj_' num2str(num_par_inj,'%1.0e')...
    '_max_time_' num2str(max_time,'%1.0e')...
    '_time_step_' num2str(step_inc,'%1.0e')...
    '_bgb_' num2str(icb.gb,'%1.0e')...
    '_gb_' num2str(ic.gb,'%1.0e')];

load(fullfile(data_save_folder,[str '.mat']))
%% 📊 Plotting the Distribution Function from Repeated Injections

% --- Figure configuration ---
FigureSize = [0 0 21 13];                  % 📐 Size of the figure in centimeters
DefaultFontSizeForFigure = 14;            % 🔤 Default font size for axes

fig = figure('Units', 'centimeters', ...
    'Position', FigureSize, ...
    'DefaultAxesFontSize', DefaultFontSizeForFigure);  % 🖼️ Create the figure

% --- Plot control parameters ---
num_of_plots = 8;                          % 🎞️ Number of snapshots to display
len = max(size(data.f_sim));              % 📏 Total number of injections (frames)
ind = [1 floor(len/num_of_plots : len/num_of_plots : len)];  % 📍 Selected indices to plot

% --- Build the constant injection steady-state distribution ---
% Each snapshot in data.f_sim represents the evolving distribution of 
% a single monoenergetic injection, recorded at time t_i during its thermalization. 
% The constant-injection DF is built as a sum over all such contributions
constant_injection_DF = zeros(size(data.f_sim(1,:)));  % 🧪 Initialize total DF

for i = 1:len
    % Normalize each individual DF before summing (area = 1)
    constant_injection_DF = constant_injection_DF + ...
        data.f_sim(i,:) / trapz(data.bins, data.f_sim(i,:));
end

% --- Plot the total constant-injection distribution ---
plot(data.bins, constant_injection_DF, 'LineWidth', 2)  % 🟢 Total distribution from all injections
hold on

% --- Plot selected individual injection snapshots ---
for i = ind
    % Normalize and plot each selected single-injection snapshot
    plot(data.bins, data.f_sim(i,:) / trapz(data.bins, data.f_sim(i,:)))
    hold on
end

% --- Configure axes ---
ax = gca;
ax.XScale = 'log';                        % 🔢 Log scale on x-axis (energy)
ax.YScale = 'log';                        % 🔢 Log scale on y-axis (DF)
ax.FontSize = 18;                         % 🔤 Tick font size

% --- Label axes with LaTeX formatting ---
xlabel('$E_k$', 'Interpreter', 'latex')   % ⚡ X-axis label: kinetic energy
ylabel('$F(E_k)$', 'Interpreter', 'latex')% 📊 Y-axis label: normalized distribution

% --- Export figure ---
exportgraphics(fig, 'Asceleration_distribution_function.png', 'Resolution', 300)  % 💾 Save as PNG

