function [data] =Single_injection_evolution_ulta_rel(N,ic,icb,max_time,time_step)
% Simulates the evolution of a population of N particles undergoing 
% ultra-relativistic collisions with a background field characterized by icb.
% Energy histograms and mean energy are recorded at fixed time intervals.

% Add base functions to path (e.g., randdir, collision_freq, etc.)
addpath('Base_func');

% Define logarithmically spaced energy bin edges for histograms, 
% spanning a wide dynamic range from 1e-16 to 1e16
edges=10.^(0:0.1:40);

%% Initialize random directions and energies
[px,py,pz,energ] = ChooseRandomRelV(ic.Energy0,N);

%% Initialize counters and parameters
col_done=0;             % Total number of collisions
step=0;                 % Next output time
time=0;                 % Total elapsed time
f_max=0;                % Maximal observed collision frequency
idx=0;                  % Output index
Energy_mean=ic.Energy0; % Initial mean energy
col_prev=0;              % Collisions counted at last output
  
% Preallocate structure to store simulation data at each output step
max_steps = ceil(max_time / time_step);
data = struct('EnergyBE1',cell(1,max_steps), ...
           'DFE1',cell(1,max_steps), ...
           'coldone',cell(1,max_steps), ...
           'coltry',cell(1,max_steps), ...
           'Energ_mean',cell(1,max_steps), ...
           'time', cell(1,max_steps));
%% Main simulation loop: stop if energy explodes or time is up
while Energy_mean<1e35 && time<max_time
    % Randomly choose a particle
    pID1=randi(N);
    p1=[px(pID1) py(pID1) pz(pID1)]; 
    e1=energ(pID1);

    % Sample random direction for the incoming background particle
    beta=icb.beta*randdir;

    % Compute current collision frequency
    f_now =collision_freq(p1,energ(pID1),beta);
    % Initialize or update maximum collision frequency
    if f_max==0
        f_max=f_now;
    end
    % Advance time using the maximal frequency (Monte Carlo approach)
    time=time+(1/N)*(1/f_max);
    f_max=max([f_max f_now]);

    % Metropolis-type acceptance test for collision
    if f_now>f_max*rand
        % Compute post-collision kinematics
        [p1,e1] = Inf_Rel_Col(p1,e1,beta);
        
        % Update particle state
        px(pID1)=p1(1); 
        py(pID1)=p1(2); 
        pz(pID1)=p1(3);
        energ(pID1)=e1; 

        % Count the collision    
        col_done=col_done+1;
    end
    % Time to record statistics
    if time>step
        step=step+time_step; 
        idx=idx+1;

        % Compute histogram of energy distribution
        [data(idx).EnergyBE1,data(idx).DFE1] = ...
            histlog(energ,edges);
        data(idx).DFE1=data(idx).DFE1./N;% Normalize

        % Store number of collisions during this interval
        col_this_time_step=col_done-col_prev;
        col_prev=col_done;

        % Store state information
        data(idx).coldone=col_done;
        data(idx).time=time;
        Energy_mean=mean(energ(1:N));
        data(idx).Energ_mean=Energy_mean;
        
        % Display progress every ~10 time units
        if mod(floor(time),10) == 0
            disp({'time step' 'ts col' 'Mean energy';...
                time col_this_time_step  Energy_mean})
        end
     end
end

% Trim unused entries
data = data(1:idx);
%% Save results
folder = 'Ultra_relativistic_evolution';
if ~exist(folder, 'dir')
    mkdir(folder);
end

str = ['Single_injection_N_' num2str(N, '%.0e') ...
       '_beta_' num2str(icb.beta, '%.0e') ...
       '_E0_' num2str(ic.Energy0, '%.0e') ...
       '_step_' num2str(time_step, '%.1f')];

save(fullfile(folder, [str '.mat']), 'data', 'ic',"icb", 'edges', 'max_time',"time_step");

end

