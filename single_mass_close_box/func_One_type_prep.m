function [Enb,data,ic]=func_One_type_prep(m,gb,N,step_inc,max_time)
addpath('Basic functions/')
%% initial condition "ic" struct
b=sqrt(gb^2/(1+gb^2));
g=1/sqrt(1-b^2);

ic.m=m;
ic.e=m*g;
ic.N=N;
ic.ke=ic.e-ic.m;
ic.edges=10.^(-16:0.1:16);
%% Prepare initial energies population
Enb=ic.e*ones(1,ic.N);
% Determine the temperature
syms x
temp=3/2*(ic.e-m);
T = double(vpasolve(3*x+ic.m*besselk(1,ic.m/x)/besselk(2,ic.m/x)==ic.e,x,temp));
clear x
%% Time considiration
t=0;
step=step_inc;  idx=1;
%% Maximal frequency and cross section dependencies
f_max=0; % start with 0 maximal frequency 
cs_p=0; % \sigma = E^-cs_p;
%% Random diraction
dir_mat = randdir_matrix(N);
ind_dir = 1;
%% data preparation
trys=0; coll_num=0; 
temp1=0; temp2=0;
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
    id1=randi(ic.N);
    id2=randi(ic.N);
    while id2==id1
        id2=randi(ic.N);
    end

    % colision prepare
    e1=Enb(id1);
    e2=Enb(id2);
    p1=sqrt((Enb(id1))^2-ic.m^2)*dir_mat(ind_dir,:);
    p2=sqrt((Enb(id2))^2-ic.m^2)*dir_mat(ind_dir+1,:);
    ind_dir=mod(ind_dir+2,N);
    if ind_dir==0
        ind_dir=1;
    end
    f_now =1e1 * collision_freq(p1,e1,p2,e2,cs_p);

    f_max = max([f_now f_max]);

    t=t+1/(ic.N/2)*(1/f_max);
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
        data(idx).kin_energy=mean(Enb-ic.m);
        [data(idx).bins_Ek,data(idx).f_sim_Ek] = histlog(Enb-ic.m,ic.edges);
        data(idx).f_sim_Ek=data(idx).f_sim_Ek/ic.N;
        data(idx).Teo_Ek=MJDF(data(idx).bins_Ek,ic.m,T);
        
        % reduced_momentum
        mu=sqrt((Enb./ic.m).^2-1); %
        [data(idx).bins_mu,data(idx).f_sim_mu] = histlog(mu,ic.edges);
        data(idx).f_sim_mu=data(idx).f_sim_mu/ic.N;
        data(idx).Teo_mu=MJmu(data(idx).bins_mu,ic.m,T);
        
        data(idx).coll_num=coll_num/N;
        data(idx).time=t;
        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
            t (data(idx).kin_energy) trys-temp1 coll_num-temp2;...
            'f_max' 'max_e' 'min_e' 'ind_dir' ;...
            f_max max(Enb-ic.m) min(Enb-ic.m) ind_dir;})
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
mkdir('saves')
save(fullfile("saves",[str '.mat']),'Enb',"data","ic");
end
