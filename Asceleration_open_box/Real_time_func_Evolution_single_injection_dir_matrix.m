function [data,str]=Real_time_func_Evolution_single_injection_dir_matrix(f_param,num_par_inj,max_time,cs_p,step_inc)
addpath('Basic functions/')
%% initial condition
ic.m=1;
ic.e=ic.m/2;
ic.edges=10.^(-16:0.1:16);
%% prepare background energies
icb.m=1e6;
icb.e=icb.m/2;
icb.N=1e6;

dir_N=1e4;
dir_mat=randdir_matrix(dir_N);
ind_dir=1;

str=['mass_' num2str(icb.m,'%1.0e')...
    '_energy_' num2str(icb.e,'%1.0e')...
    '_particle_number_' num2str(icb.N,'%1.0e')];
folder='saves_one_type';
datab=load(fullfile(folder,[str '.mat']));
En=datab.Energy_array;
%% first injection
ic.m=1;
ic.e=ic.m/2;
ic.N=num_par_inj;

str=['mass_' num2str(ic.m,'%1.0e')...
     '_energy_' num2str(ic.e,'%1.0e')...
     '_particle_number_' num2str(ic.N,'%1.0e')];

folder='saves_one_type';
Array=load(fullfile(folder,[str '.mat']));

Energy_low_mass=Array.Energy_array;
num_of_par=num_par_inj;
N=num_of_par;
T = 2/3*(icb.e);
%% time
t=0;
% step_inc=0.1;
step=0; idx=1;
%% maximal frequancy and cross section dependencies
f_max=0; dt=0;
% cs_p=0;
%% data preparation
trys=0; coll_num=0; temp1=0; temp2=0;
data = struct('f_sim', cell(1, max_time),...
    'bins', cell(1, max_time),...
    'kin_energy', cell(1, max_time));
%% MB func
MBDF=@(e,T) 2/(sqrt(pi*T^3)).*exp(-e/T).*sqrt(e);
%% Collision loop
tic
while t<max_time
    if t>=step
        data(idx).kin_energy=mean(Energy_low_mass);
         [data(idx).bins,data(idx).f_sim] = histlog(Energy_low_mass,ic.edges);
        data(idx).f_sim=data(idx).f_sim/trapz(data(idx).bins,data(idx).f_sim);

%         data(idx).Teo=MBDF(data(idx).bins,T);
        data(idx).time=t;

        % new diractions
        dir_mat=randdir_matrix(dir_N);
% 
%         disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
%             t data(idx).kin_energy trys-temp1 coll_num-temp2; ...
%             'f_max' 'num_of_par' 'max_energy' 'dt :)'; ...
%             f_max num_of_par max(Energy_low_mass) dt})
%         temp1=trys;
%         temp2=coll_num;
        step=step+step_inc;
        idx=idx+1;
        N=num_of_par;
    end
    trys=trys+1;

    % index selection
    id1=randi(num_of_par);
    id2=randi(icb.N);

    % colision prepare
    v1=sqrt(2*Energy_low_mass(id1)/ic.m)*dir_mat(ind_dir,:);
    v2=sqrt(2*En(id2)/icb.m)*dir_mat(ind_dir+1,:);
    ind_dir=mod(ind_dir+2,dir_N);
    if ind_dir==0
        ind_dir=1;
    end
    f_now =f_param*frequancy(v1,v2,cs_p);
    f_max = max([f_now f_max]);
    dt=1/(N)*(1/f_max);
    t=t+dt;
    if f_now > f_max*rand
        coll_num=coll_num+1;
        n0=dir_mat(ind_dir,:);
        ind_dir=mod(ind_dir+1,dir_N);
        if ind_dir==0
            ind_dir=1;
        end
        [Energy_low_mass(id1),~] = collision(v1,v2,ic.m,icb.m,n0);
    end
    
end
toc
folder='saves_v4';
mkdir(folder);
str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
    '_step_inc_' num2str(step_inc,'%.2f')];

% save(fullfile(folder,[str '_mat_dir.mat']));
% anima_s(data,str,folder);
end