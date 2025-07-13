function [Enb,data,icb]=func_One_type_prep(m,gb,N,step_inc)
addpath('Basic functions/')
%% initial condition

% gb=1e-2;
b=sqrt(gb^2/(1+gb^2));
g=1/sqrt(1-b^2);

icb.m=m;
icb.e=m*g;
icb.N=N;
icb.ke=icb.e-icb.m;
icb.edges=10.^(-16:0.1:16);
%% prepare initial energies
Enb=icb.e*ones(1,icb.N);
syms x
temp=3/2*(icb.e-m);
T = double(vpasolve(3*x+icb.m*besselk(1,icb.m/x)/besselk(2,icb.m/x)==icb.e,x,temp));
clear x
%% time
t=0;
step=step_inc;  idx=1;
max_time=6;
%% maximal frequancy and cross section dependencies
f_max=0;
cs_p=0;
%% Random diraction
dir_mat = randdir_matrix(N);
ind_dir = 1;
%% data preparation
trys=0; coll_num=0; temp1=0; temp2=0;
data = struct('f_sim', cell(1, max_time),...
    'bins', cell(1, max_time),...
    'kin_energy', cell(1, max_time));
%% MB func
MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
    .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);
MJmu=@(x,n,T) (m/T)*(1/besselk(2,m/T)).*exp(-(m/T)*sqrt(1+x.^2)).*x.^2;
%% Collision loop
tic
while t<max_time
    trys=trys+1;
    % index selection
    id1=randi(icb.N);
    id2=randi(icb.N);
    while id2==id1
        id2=randi(icb.N);
    end

    % colision prepare
    e1=Enb(id1);
    e2=Enb(id2);
    p1=sqrt((Enb(id1))^2-icb.m^2)*dir_mat(ind_dir,:);
    p2=sqrt((Enb(id2))^2-icb.m^2)*dir_mat(ind_dir+1,:);
    ind_dir=mod(ind_dir+2,N);
    if ind_dir==0
        ind_dir=1;
    end
    f_now =1e1 * collision_freq(p1,e1,p2,e2,cs_p);

    f_max = max([f_now f_max]);

    t=t+1/(icb.N/2)*(1/f_max);
    if f_now > f_max*rand
        coll_num=coll_num+1;
        n0=dir_mat(ind_dir+1,:);
        ind_dir=mod(ind_dir+1,N);
        if ind_dir==0
            ind_dir=1;
        end

        [Enb(id1), ~,Enb(id2), ~] = RelCol(e1, p1,e2, p2, n0);
    end
    if t>step
        data(idx).kin_energy=mean(Enb-icb.m);
        [data(idx).bins_Ek,data(idx).f_sim_Ek] = histlog(Enb-icb.m,icb.edges);
        data(idx).f_sim_Ek=data(idx).f_sim_Ek/icb.N;
        data(idx).Teo_Ek=MJDF(data(idx).bins_Ek,icb.m,T);
        
        mu=sqrt((Enb./icb.m).^2-1);
        [data(idx).bins_mu,data(idx).f_sim_mu] = histlog(mu,icb.edges);
        data(idx).f_sim_mu=data(idx).f_sim_mu/icb.N;
        data(idx).Teo_mu=MJmu(data(idx).bins_mu,icb.m,T);
        
        data(idx).coll_num=coll_num/N;
        data(idx).time=t;
        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
            t (data(idx).kin_energy) trys-temp1 coll_num-temp2;...
            'f_max' 'max_e' 'min_e' 'ind_dir' ;...
            f_max max(Enb-icb.m) min(Enb-icb.m) ind_dir;})
        temp1=trys;
        temp2=coll_num;
        step=step+step_inc;
        idx=idx+1;
        %% Random diraction
        dir_mat = randdir_matrix(N);
        ind_dir = 1;
    end
end
toc
%%
str=['One_type_mass_' num2str(m,'%1.0e')...
    '_initial_gb_' num2str(gb,'%1.0e')...
    '_particle_number_' num2str(N,'%1.0e')];
% anima_s(data,str,'Videos1');
save(fullfile("saves",[str '.mat']),'Enb',"data","icb");
end