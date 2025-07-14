clear
clc
addpath('Basic functions/')
%% parameters
max_time=1e3;
for cs_p=[0]
    for step_inj_inc=[0.1 1 5 10 20 100] 
        for energy_inc=[1]
            for nop_energy_inc=100
                %% initial condition
                ic.m=1;
                ic.N=1e6;
                ic.edges=10.^(-16:0.1:16);
                ic.gb=0.05;
                ic.g=sqrt(ic.gb^2+1)
                ic.e=ic.m*ic.g;
                ic.ek=ic.e-ic.m;
                %% Prepare initial energies
                % sum(En)
                % En=ic.e*ones(1,ic.N);
                % load('s.mat','En')
                %                 En=func_One_type_prep(ic.m,ic.gb,ic.N);

                load(fullfile('Background_population',['MJ_Mass_' num2str(ic.m,'%1.0e') '_gb_' num2str(ic.gb,'%1.0e')...
                    '_partcle_number_' num2str(ic.N,'%1.0e')]))

                syms x
                T = double(vpasolve(3*x+ic.m*besselk(1,ic.m/x)/besselk(2,ic.m/x)==ic.e,x));
                clear x
                %% time
                t=0; f_param=1;
                step_inc=1; step=step_inc; idx=1;
                step_inj_inc_one_particle=step_inj_inc/nop_energy_inc;
                step_inj=step_inj_inc_one_particle;

                %% maximal frequancy and cross section dependencies
                f_max=0;

                %% data preparation
                trys=0; coll_num=0; temp1=0; temp2=0;

                data = struct('f_sim', cell(1, max_time),...
                    'bins', cell(1, max_time),...
                    'kin_energy', cell(1, max_time));
                %% diraction matrix
                N_dir=1e6;
                Mat_dir=randdir_mat(N_dir);
                dir_ind=1;
                %% MB func
                MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
                    .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);
                %% Collision loop
                tic
                while t<max_time
                    trys=trys+1;
                    % index selection
                    id1=randi(ic.N);
                    id2=randi(ic.N);
                    while id2==id1
                        id2=randi(ic.N);
                    end
                    % colision prepare
                    e1=En(id1);
                    e2=En(id2);
                    p1=sqrt(En(id1)^2-ic.m^2)*Mat_dir(dir_ind,:);
                    if dir_ind==N_dir
                        dir_ind=1;
                    else
                        dir_ind=dir_ind+1;
                    end
                    p2=sqrt(En(id2)^2-ic.m^2)*Mat_dir(dir_ind,:);

                    if dir_ind==N_dir
                        dir_ind=1;
                    else
                        dir_ind=dir_ind+1;
                    end

                    f_now= f_param*collision_freq(p1,e1,p2,e2,cs_p);
                    f_max = max([f_now f_max]);

                    t=t+1/(ic.N/2)*(1/f_max);
                    if f_now > f_max*rand
                        coll_num=coll_num+1;
                        n0=Mat_dir(dir_ind,:);
                        if dir_ind==N_dir
                            dir_ind=1;
                        else
                            dir_ind=dir_ind+1;
                        end
                        [En(id1), p1_after,En(id2), p2_after] =...
                            RelCol(e1, p1,e2, p2,n0);
                    end
                    if t>step
                        %         toc
                        %         tic
                        data(idx).kin_energy=mean(En-ic.m);
                        [data(idx).bins,data(idx).f_sim] = histlog(En-ic.m,ic.edges);
                        data(idx).f_sim=data(idx).f_sim/ic.N;
                        data(idx).Teo=MJDF(data(idx).bins,ic.m,T);
                        data(idx).time=t;
                        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
                            t data(idx).kin_energy trys-temp1 coll_num-temp2; ...
                            'f_max' 'Total_energy' 'np_above_ke' 'max energy'; ...
                            f_max sum(En-ic.m) sum((En-ic.m)>data(idx).kin_energy) max(En-ic.m)})
                        temp1=trys;
                        temp2=coll_num;
                        step=step+step_inc;
                        %% time update
                        idx=idx+1;
                           %% diraction matrix
                        N_dir=1e6;
                        Mat_dir=randdir_mat(N_dir);
                        dir_ind=1;
                    end
                    if t>step_inj
                        step_inj=step_inj+step_inj_inc_one_particle;
                     %%  inject energy to one particle
                        inc_ind=randi(ic.N,1,1); % choose the one particle for assceleration
                        for idr=inc_ind
                            % energy inc
                            En(idr)=En(idr)+energy_inc;
                            % random particle cooling

                            if prod(En>=0)==0
                                error('negative energy')
                            end

                            Enk=En-ic.m;
                            Enk=(ic.N*(ic.e-ic.m))/(ic.N*(ic.e-ic.m)+energy_inc)*Enk;
                            En=Enk+ic.m;
                        end
                    end
                end
                toc
                %%
                folder = 'saves_one_by_one_inj_v1';
                mkdir(folder);
                str=['N_' num2str(ic.N,'%1.0e')...
                    '_cs_p_'  num2str(cs_p)...
                    '_energy_inc_' num2str(energy_inc,'%1.0e')...
                    '_max_time_' num2str(max_time,'%1.0e')...
                    '_nop_par_each_ts_'  num2str(nop_energy_inc,'%1.0e')...
                    '_step_inj_inc_' num2str(step_inj_inc,'%1.0e') ];
                save([fullfile(folder,str) '.mat'],'-mat')
                % anima_s(data,str,folder)
            end
        end
    end
end
%%
% cftool(log10(data(end).bins),log10(data(end).f_sim))