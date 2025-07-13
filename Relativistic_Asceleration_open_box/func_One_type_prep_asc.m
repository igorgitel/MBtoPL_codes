function [data En]=func_One_type_prep_asc(f_param,num_par_inj,max_time)
addpath('Basic functions/')
%% background load 
m=1e9;
bgb=1e-1;
N=1e5;

b=sqrt(bgb^2/(1+bgb^2));
g=1/sqrt(1-b^2);
eb=m*g;

str=['One_type_mass_' num2str(m,'%1.0e')...
    '_initial_gb_' num2str(bgb,'%1.0e')...
    '_particle_number_' num2str(N,'%1.0e')];
%  [Enb,~,icb]=func_One_type_prep(m,bgb,N)
 load(fullfile('saves',[str '.mat']),"Enb","icb")
BKE=eb-m;
%% initial condition
gb=1e-4;
b=sqrt(gb^2/(1+gb^2));
g=1/sqrt(1-b^2);

ic.m=1;
ic.e=g*ic.m;


KE=ic.e-ic.m;
ic.edges=10.^(-16:0.1:16);
disp({'Back_kin_energ' 'par_kin_enrg';
     BKE KE})
clear m N
%% prepare initial energies
En=ic.e*ones(1,num_par_inj);
num_par=num_par_inj;
syms x
T_temp=3/2*(icb.e-icb.m);
T = double(vpasolve(3*x+icb.m*besselk(1,icb.m/x)/besselk(2,icb.m/x)==icb.e,x,T_temp));
clear x
%% time
t=0;
step=1; step_inc=1; idx=1;
%% maximal frequancy and cross section dependencies
f_max=0;
cs_p=0;
%% Random diraction 
N_dir=1e5;
dir_mat = randdir_matrix(N_dir);
ind_dir = 1;
%% data preparation
trys=0; coll_num=0; temp1=0; temp2=0;
data = struct('f_sim', cell(1, max_time),...
    'bins', cell(1, max_time),...
    'kin_energy', cell(1, max_time));
%% MB func
MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
        .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);
%% Collision loop
tic
while t<max_time
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

    t=t+1/(num_par)*(1/f_max);
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
        data(idx).kin_energy=mean(En-ic.m);
        [data(idx).bins,data(idx).f_sim] = histlog(En-ic.m,ic.edges);
        data(idx).f_sim=data(idx).f_sim/num_par;
        data(idx).Teo=MJDF(data(idx).bins,icb.m,T);
        data(idx).time=t;
        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
            t (data(idx).kin_energy) trys-temp1 coll_num-temp2;...
            'f_max' 'BKE' 'KE' 'num_par' ;...
            f_max BKE KE num_par;})
        temp1=trys;
        temp2=coll_num;
        step=step+step_inc;
        idx=idx+1;
        %% Random diraction 
        dir_mat = randdir_matrix(N_dir);
        ind_dir = 1;
        %% inject new particles
        temp=num_par+num_par_inj;
        En((num_par+1):temp)=ic.e*ones(1,num_par_inj);
        num_par=temp;
    end
end
toc
% data=data(1:5860)
%%
 str=['Constant_injection_f_param_' num2str(f_param,'%1.0e')...
     '_num_par_inj_' num2str(num_par_inj,'%1.0e')...
     '_max_time_' num2str(max_time,'%1.0e')...
     '_Bgb_' num2str(bgb,'%1.0e')...
     '_gb_' num2str(gb,'%1.0e')];
anima_s(data,str)
save(fullfile("saves",[str '.mat']),'data','En',"ic");
end