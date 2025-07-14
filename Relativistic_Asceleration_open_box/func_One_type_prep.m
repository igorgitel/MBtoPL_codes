function [Enb,data,icb]=func_One_type_prep(m,gb,N,folder)
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
step_inc=1; step=step_inc;  idx=1;
max_time=1e4;
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
        [data(idx).bins,data(idx).f_sim] = histlog(Enb-icb.m,icb.edges);
        data(idx).f_sim=data(idx).f_sim/icb.N;

        data(idx).Teo=MJDF(data(idx).bins,icb.m,T);
        data(idx).time=t;
        data(idx).ssd = sum((data(idx).f_sim - data(idx).Teo).^2); %sum of squared deviations
        data(idx).mse = mean((data(idx).f_sim - data(idx).Teo).^2);  % Mean Squared Error
        disp({'ts' 'mean energy' 'trys_c' 'coll_num_c'; ...
            t (data(idx).kin_energy) trys-temp1 coll_num-temp2;...
            'f_max' 'max_e' 'min_e' 'mse' ;...
            f_max max(Enb-icb.m) min(Enb-icb.m) data(idx).mse;})
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
mkdir(folder)
str=['One_type_mass_' num2str(m,'%1.0e')...
    '_initial_gb_' num2str(gb,'%1.0e')...
    '_particle_number_' num2str(N,'%1.0e')];
icb.str=str;
% anima_s(data,str,folder);

idx = numel(data);  % or simply: idx = end;

idx = numel(data);  % last index

% Extract values
x = data(idx).bins;
y = data(idx).f_sim;
y_teo = data(idx).Teo;

% Compute dynamic limits
x_min = min(x(y > 0));
x_max = max(x(y > 0));
y_min = min(y(y > 0));
y_max = max(y(y > 0));

% Expand by one order of magnitude
xlim_vals = [10^floor(log10(x_min)), 10^ceil(log10(x_max))];
ylim_vals = [10^floor(log10(y_min)), 10^ceil(log10(y_max))];

% Plot
figure;
loglog(x, y, 'o-', 'LineWidth', 1.5, 'DisplayName', 'Simulated');
hold on;
loglog(x, y_teo, '-', 'LineWidth', 1.5, 'DisplayName', 'Theoretical');
hold off;

xlim(xlim_vals);
ylim(ylim_vals);

xlabel('Energy');
ylabel('f(E)');
legend('Location', 'northeast');
grid on;
title(['Comparison at t = ', num2str(data(idx).time)]);

save(fullfile(folder,[str '.mat']),'Enb',"icb");
end