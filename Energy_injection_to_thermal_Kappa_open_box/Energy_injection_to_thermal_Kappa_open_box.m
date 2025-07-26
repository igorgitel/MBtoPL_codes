clear
clc
addpath('Basic functions/')

%% parameters
max_time=1e3; % total simulation time

for cs_p=[0] % cross-section power-law exponent (sigma = E^(-cs_p))
    for step_inj_inc=[10 100 1000]  % time interval between energy injections
        for energy_inc=[1] % amount of energy injected per event
            for nop_energy_inc = [100] % number of injections per time unit

                %% initial condition
                ic.m=1; % particle mass
                ic.N=1e6; % number of particles
                ic.edges=10.^(-16:0.1:16); % histogram bin edges (log-spaced)
                ic.gb=0.05; % initial beta
                ic.g=sqrt(ic.gb^2+1); % corresponding gamma factor
                ic.e=ic.m*ic.g; % total initial energy
                ic.ek=ic.e-ic.m; % initial kinetic energy

                %% Prepare initial energies (load or simulate)
                str=['One_type_mass_' num2str(ic.m,'%1.0e')...
                    '_initial_gb_' num2str(ic.gb,'%1.0e')...
                    '_particle_number_' num2str(ic.N,'%1.0e')];
                background_file_path=fullfile('saves_one_type',[str '.mat']);

                if exist(background_file_path, 'file') == 2
                    load(background_file_path)
                else
                    step_inc_back = 1;
                    max_time_back = 10;
                    [Energy_array,data,~] = func_One_type_prep(ic.m,ic.gb,ic.N,step_inc_back,max_time_back);
                end

                %% calculate temperature from inverse Maxwell-Juttner equation
                syms x
                T = double(vpasolve(3*x+ic.m*besselk(1,ic.m/x)/besselk(2,ic.m/x)==ic.e,x));
                clear x

                %% Optional: visualize initial distribution
                check_initial_condition = false;
                if check_initial_condition
                    last_idx = find(~cellfun(@isempty, {data.kin_energy}), 1, 'last');
                    Ek = data(last_idx).bins_Ek;
                    f_sim = data(last_idx).f_sim_Ek;
                    f_teo = data(last_idx).Teo_Ek;

                    figure;
                    loglog(Ek, f_sim, '-o'); hold on;
                    loglog(Ek, f_teo, '-');
                    xlabel('$E_k$', 'Interpreter', 'latex');
                    ylabel('Distribution function', 'Interpreter', 'latex');
                    legend({'$f_{\mathrm{sim}}$', '$f_{\mathrm{theory}}$'}, 'Interpreter', 'latex', 'Location', 'best');
                    ylim([1e-4,1e4])
                    xlim([1e-6,1e-1])
                    grid on;
                    title('Simulated vs. Theoretical Energy Distribution');
                end

                %% time initialization
                t=0; f_param=1;
                step_inc=1; step=step_inc; idx=1;
                step_inj_inc_one_particle=step_inj_inc/nop_energy_inc;
                step_inj=step_inj_inc_one_particle;

                %% maximum frequency tracking
                f_max=0;

                %% preallocate data structure
                trys=0; coll_num=0; temp1=0; temp2=0;
                data = struct('f_sim', cell(1, max_time),...
                    'bins', cell(1, max_time),...
                    'kin_energy', cell(1, max_time));

                %% random direction matrix for collisions
                N_dir=1e6;
                Mat_dir=randdir_matrix(N_dir);
                dir_ind=1;

                %% Maxwell-Juttner distribution function
                MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
                    .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);

                %% Main collision loop
                tic
                while t<max_time
                    trys=trys+1;

                    % randomly select two different particles
                    id1=randi(ic.N);
                    id2=randi(ic.N);
                    while id2==id1
                        id2=randi(ic.N);
                    end

                    % get energy and momentum vectors
                    e1=Energy_array(id1);
                    e2=Energy_array(id2);
                    p1=sqrt(Energy_array(id1)^2-ic.m^2)*Mat_dir(dir_ind,:);
                    dir_ind = mod(dir_ind,N_dir) + 1;
                    p2=sqrt(Energy_array(id2)^2-ic.m^2)*Mat_dir(dir_ind,:);
                    dir_ind = mod(dir_ind,N_dir) + 1;

                    % compute collision frequency
                    f_now= f_param*collision_freq(p1,e1,p2,e2,cs_p);
                    f_max = max([f_now f_max]);

                    % time step update
                    t=t+1/(ic.N/2)*(1/f_max);

                    % collision acceptance
                    if f_now > f_max*rand
                        coll_num=coll_num+1;
                        n0=Mat_dir(dir_ind,:);
                        dir_ind = mod(dir_ind,N_dir) + 1;
                        [Energy_array(id1), ~, Energy_array(id2), ~] =...
                            RelCol(e1, p1,e2, p2,n0);
                    end

                    % store state at intervals
                    if t>step
                        data(idx).kin_energy=mean(Energy_array-ic.m);
                        [data(idx).bins,data(idx).f_sim] = histlog(Energy_array-ic.m,ic.edges);
                        data(idx).f_sim=data(idx).f_sim/ic.N;
                        data(idx).Teo=MJDF(data(idx).bins,ic.m,T);
                        data(idx).time=t;
                        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
                            t data(idx).kin_energy trys-temp1 coll_num-temp2; ...
                            'f_max' 'Total_energy' 'np_above_ke' 'max energy'; ...
                            f_max sum(Energy_array-ic.m) sum((Energy_array-ic.m)>data(idx).kin_energy) max(Energy_array-ic.m)})
                        temp1=trys;
                        temp2=coll_num;
                        step=step+step_inc;

                        % regenerate direction matrix
                        N_dir=1e6;
                        Mat_dir=randdir_matrix(N_dir);
                        dir_ind=1;
                    end

                    % inject energy to one particle
                    if t>step_inj
                        step_inj=step_inj+step_inj_inc_one_particle;
                        inc_ind=randi(ic.N,1,1); % random particle
                        for idr=inc_ind
                            Energy_array(idr)=Energy_array(idr)+energy_inc;

                            if prod(Energy_array>=0)==0
                                error('negative energy')
                            end

                            % rescale all other particles (cooling)
                            Energy_arrayk=Energy_array-ic.m;
                            Energy_arrayk=(ic.N*(ic.e-ic.m))/(ic.N*(ic.e-ic.m)+energy_inc)*Energy_arrayk;
                            Energy_array=Energy_arrayk+ic.m;
                        end
                    end
                end
                toc

                %% Save results and generate animation
                folder = 'saves_one_by_one_inj_v1';
                mkdir(folder);
                str=['N_' num2str(ic.N,'%1.0e')...
                    '_cs_p_'  num2str(cs_p)...
                    '_energy_inc_' num2str(energy_inc,'%1.0e')...
                    '_max_time_' num2str(max_time,'%1.0e')...
                    '_nop_par_each_ts_'  num2str(nop_energy_inc,'%1.0e')...
                    '_step_inj_inc_' num2str(step_inj_inc,'%1.0e') ];
                save([fullfile(folder,str) '.mat']);
                anima_s(data,str,folder)
            end
        end
    end
end
