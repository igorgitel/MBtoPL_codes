function [data,En]=func_asc_single_injection_Mat_TS( ...
    ic,icb,f_param,num_par_inj,max_time,step_inc,data_save_folder)
addpath('Basic functions/')
%% background load

str=['One_type_mass_' num2str(icb.m,'%1.0e')...
    '_initial_gb_' num2str(icb.gb,'%1.0e')...
    '_particle_number_' num2str(icb.N,'%1.0e')];

load(fullfile('saves_one_type',[str '.mat']),"Energy_array")
Enb=Energy_array; % assign for background 
% add properties
icb.b=sqrt(icb.gb^2/(1+icb.gb^2));
icb.g=1/sqrt(1-icb.b^2);
icb.e=icb.m*icb.g;
icb.ke = icb.e - icb.m;
clear Energy_array

%edges definition
edges=10.^(-16:0.1:16);
%% The background dist

MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
    .*(sqrt((x+m).^2-m.^2)).*exp(-x/T); % MJ func
% Temp calc
syms x
temp=3/2*(icb.e-icb.m);
T = double(vpasolve(3*x+...
    icb.m*besselk(1,icb.m/x)/besselk(2,icb.m/x)==icb.e,x,temp));

%% initial condition
ic.b=sqrt(ic.gb^2/(1+ic.gb^2));
ic.g=1/sqrt(1-ic.b^2);

ic.m=1;
ic.e=ic.g*ic.m;
ic.ke=ic.e-ic.m;

disp({'Back_kin_energ' 'par_kin_enrg';
    icb.ke ic.ke})
%% prepare initial energies

En=ic.e*ones(1,num_par_inj);
num_par=num_par_inj;

%% time
t=0;
step=step_inc;  idx=1;
show_step_inc=1; show_step=show_step_inc;

%% maximal frequancy and cross section dependencies
f_max=0;
cs_p=0;

%% Random diraction
N_dir=1e5;
dir_mat = randdir_matrix(N_dir);
ind_dir = 1;
%% data preparation
trys=0; coll_num=0; temp1=0; temp2=0;
KE_condition=ic.ke;
% data = struct('f_sim',...
    % 'bins',...
    % 'kin_energy');
data.f_sim=zeros(max_time/step_inc,(length(edges)-1));
data.kin_energy=zeros(max_time/step_inc,1);
data.time=zeros(max_time/step_inc,1);
%% Collision loop

while t < max_time && KE_condition <= icb.ke
  
    trys=trys+1;
    % index selection
    id1=randi(num_par);
    id2=randi(icb.N);

    % colision prepare
    e1=En(id1);
    e2=Enb(id2);
    p1=sqrt((En(id1))^2-ic.m^2)*dir_mat(ind_dir,:);
    p2=sqrt((Enb(id2))^2-icb.m^2)*dir_mat(ind_dir+1,:);
    ind_dir=mod(ind_dir+2,N_dir);
    if ind_dir==0
        ind_dir=1;
    end
    f_now =f_param * collision_freq(p1,e1,p2,e2,cs_p);

    f_max = max([f_now f_max]);
    dt=1/(num_par)*(1/f_max);
    t=t+dt;
    if f_now > f_max*rand
        coll_num=coll_num+1;
        n0=dir_mat(ind_dir+1,:);
        ind_dir=mod(ind_dir+1,N_dir);
        if ind_dir==0
            ind_dir=1;
        end
        [En(id1), ~,~,~] = RelCol(e1, p1,e2, p2,n0);
    end
    if t>step
        data.kin_energy(idx)=mean(En-ic.m);
        KE_condition=data.kin_energy(idx);
        data.std(idx)=std(En-ic.m);
        [data.bins,data.f_sim(idx,:),data.f_sim_50(idx,:),data.f_sim_100(idx,:)]...
            = low_err_histlog(En-ic.m,edges);

        data.time(idx)=t;
        step=step+step_inc;
        
        data.coll_num_pp(idx)=coll_num/num_par;
        
        idx=idx+1;
        %% Random diraction
        dir_mat = randdir_matrix(N_dir);
        ind_dir = 1;        
        if t>show_step
                disp({'ts' 'mean energy' 'trys_PP_PerTS' 'coll_PP_PerTS'; ...
                num2str(t,'%.2f') (data.kin_energy(idx-1)) (trys-temp1)/num_par (coll_num-temp2)/num_par;...
                'f_max'  'KE' 'BKE' 'dt' ;...
                f_max ic.ke icb.ke dt;})
                show_step=show_step+show_step_inc;
        end
        temp1=trys;
        temp2=coll_num;
        %         %% inject new particles
        %         temp=num_par+num_par_inj;
        %         En((num_par+1):temp)=ic.e*ones(1,num_par_inj);
        %         num_par=temp;
    end


end
% add theoretical data 
data.Teo=MJDF(data.bins,ic.m,T);
% adjust the final structure if needed 
final_idx = idx - 1;
data.kin_energy = data.kin_energy(1:final_idx);
data.std = data.std(1:final_idx);
data.f_sim =  data.f_sim(1:final_idx,:);
data.f_sim_50 =  data.f_sim_50(1:final_idx,:);
data.f_sim_100 =  data.f_sim_100(1:final_idx,:);
data.time = data.time(1:final_idx);
data.coll_num_pp = data.coll_num_pp(1:final_idx);

%%
mkdir(data_save_folder)

str=['mass_Background_' num2str(icb.m,'%1.0e')...
    '_Single_injection_f_param_' num2str(f_param,'%1.0e')...
    '_num_par_inj_' num2str(num_par_inj,'%1.0e')...
    '_max_time_' num2str(max_time,'%1.0e')...
    '_time_step_' num2str(step_inc,'%1.0e')...
    '_bgb_' num2str(icb.gb,'%1.0e')...
    '_gb_' num2str(ic.gb,'%1.0e')];

save(fullfile(data_save_folder,[str '.mat']),'-v7.3', ...
    'ic','icb', ...
    'data','En', ...
    'max_time','step_inc', ...
    'cs_p');
% anima_s(data,str,folder);
end