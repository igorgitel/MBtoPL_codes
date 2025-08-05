function data = Energy_injection_no_cooling(ic, ...
                                           cs_p, ...
                                           energy_inc, ...
                                           num_of_inj_particles,...
                                           time_between_injections, ...
                                           max_time, ...
                                           save_folder, ...
                                           is_save,parallel_computation)
    % Function to simulate relativistic particle collisions with energy injection
    % into a single particle at regular intervals.
    % 
    % INPUTS:
    %   ic                          - struct with fields: m, N, gb
    %   cs_p                        - cross-section power-law exponent (scalar)
    %   energy_inc                  - energy injected per step (scalar)
    %   time_between_injections         - interval between injections (scalar)
    %   max_time                    - total simulation time (scalar)
    %   save_folder                 - folder to save results (string)
    %
    % OUTPUT:
    %   data                        - struct containing simulation history
    
%% Compute heating rates
Heat_rate = energy_inc * num_of_inj_particles / time_between_injections;  % Total heating rate
h = Heat_rate / ic.N;  % Specific heating rate per particle

if ~parallel_computation
    disp(['Total Heating rate = ' num2str(Heat_rate) ' energy units per time unit'])
    disp(['Heating rate h = ' num2str(h,'%1.2e') ' per particle per time unit'])
end

    if nargin < 6
        save_folder = '';
    end
    edges=10.^(-16:0.1:16);
    % Derived initial parameters
    ic.g = sqrt(ic.gb^2 + 1);
    ic.e = ic.m * ic.g;
    ic.ek = ic.e - ic.m;
        
    % Load background energy distribution
    load(fullfile('Background_population', ...
        ['MJ_Mass_' num2str(ic.m, '%1.0e') '_gb_' num2str(ic.gb, '%1.0e') ...
        '_partcle_number_' num2str(ic.N, '%1.0e')]), 'En');

    % Solve for temperature of MJ distribution
    % syms x
    % T = double(vpasolve(3*x + ic.m*besselk(1,ic.m/x)/besselk(2,ic.m/x) == ic.e, x));
    % clear x

    % Initialization
    t = 0; f_param = 1; 
    step_inc = 1; step = step_inc; idx = 1;
    step_inj = time_between_injections;

    f_max = 0; trys = 0; coll_num = 0; temp1 = 0; temp2 = 0;

    data = struct('f_sim', cell(1, max_time), ...
                  'bins', cell(1, max_time), ...
                  'kin_energy', cell(1, max_time), ...
                  'Teo', cell(1, max_time), ...
                  'time', cell(1, max_time));
    % Direction matrix initialization
    N_dir = 1e6;
    Mat_dir = randdir_mat(N_dir);
    dir_ind = 1;

    % MJDF = @(x,m,T) 1./(m.^2.*T.*besselk(2,m./T)).*(x+m).* ...
    %                sqrt((x+m).^2 - m.^2).*exp(-x./T);
    computation_time=tic;
    while t < max_time
       
        trys = trys + 1;

        % Pick two distinct particles
        id1 = randi(ic.N);
        id2 = randi(ic.N);
        while id2 == id1
            id2 = randi(ic.N);
        end

        e1 = En(id1);
        e2 = En(id2);
        p1 = sqrt(e1^2 - ic.m^2) * Mat_dir(dir_ind,:); dir_ind = mod(dir_ind, N_dir) + 1;
        p2 = sqrt(e2^2 - ic.m^2) * Mat_dir(dir_ind,:); dir_ind = mod(dir_ind, N_dir) + 1;

        f_now = f_param * collision_freq(p1,e1,p2,e2,cs_p);
        f_max = max([f_now, f_max]);

        t = t + 1/(ic.N/2) * (1/f_max);
        if f_now > f_max * rand
            coll_num = coll_num + 1;
            n0 = Mat_dir(dir_ind,:); dir_ind = mod(dir_ind, N_dir) + 1;
            [En(id1), ~, En(id2), ~] = RelCol(e1, p1, e2, p2, n0);
        end

        if t > step
            elapsed_time = toc(computation_time);
            computation_time=tic ; 
            data(idx).kin_energy = mean(En - ic.m);
            [data(idx).bins, data(idx).f_sim] = histlog(En - ic.m, edges);
            data(idx).f_sim = data(idx).f_sim / ic.N;
            % e_param = ic.e + h*t;
            % syms x
            % T = double(vpasolve(3*x + ic.m*besselk(1,ic.m/x)/besselk(2,ic.m/x) == e_param, x));
            % clear x
            % data(idx).Teo = MJDF(data(idx).bins, ic.m, T);
            data(idx).time = t;
            if ~parallel_computation
                disp({'ts','mean energy','trys','collisions'; ...
                      t, data(idx).kin_energy, trys - temp1, coll_num - temp2; ...
                      'f_max','Total_Ek','>mean Ek','Comp_time'; ...
                      f_max, sum(En - ic.m), sum((En - ic.m) > data(idx).kin_energy), elapsed_time})
            end
            temp1 = trys; temp2 = coll_num;
            step = step + step_inc;
            Mat_dir = randdir_mat(N_dir);
            dir_ind = 1;
            idx = idx + 1;
        end

        if t > step_inj
            step_inj = step_inj + time_between_injections;
            inc_ind = randi(ic.N,[1 num_of_inj_particles]);
            En(inc_ind) = En(inc_ind) + energy_inc;
            if any(En < ic.m)
                error('Negative energy encountered')
            end

            % Enk = En - ic.m;
            % Enk = (ic.N * ic.ek) / (ic.N * ic.ek + num_of_inj_particles*energy_inc) * Enk;
            % En = Enk + ic.m;
        end
    end
    
   if is_save
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        str =  ['N_' num2str(ic.N,'%1.0e') ...
               '_e_' num2str(energy_inc,'%1.0e') ...
               '_dt_' num2str(time_between_injections,'%1.0e') ...
               '_n_' num2str(num_of_inj_particles,'%1.0e') ...
               '_H_' num2str(Heat_rate,'%1.1e')...
               '_max_time_' num2str(max_time,'%1.0e') ...
               '_cs_p_' num2str(cs_p)];

        save(fullfile(save_folder, [str '.mat']), 'data', 'ic', 'cs_p', ...
            'energy_inc', 'time_between_injections', ...
            'max_time', ...
            'num_of_inj_particles', ...
            'Heat_rate');
    end
end
