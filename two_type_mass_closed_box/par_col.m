function [data, ic1,ic2] = par_col(ic1,ic2,N,max_time,fparam)
addpath('Basic functions/')
%% initial condition

ic1.b=sqrt(ic1.gb^2/(1+ic1.gb^2));
ic1.g=1/sqrt(1-ic1.b^2);
ic1.e=ic1.m*ic1.g;
ic1.ke=ic1.e-ic1.m;


ic2.b=sqrt(ic2.gb^2/(1+ic2.gb^2));
ic2.g=1/sqrt(1-ic2.b^2);
ic2.e=ic2.m*ic2.g;
ic2.ke=ic2.e-ic2.m;

%% prepare initial energies
Enb1=ic1.e*ones(1,N);
Enb2=ic2.e*ones(1,N);
Enb=[Enb1 Enb2];

%% Temperature calculation 
syms x
temp=2/3*(ic1.ke+ic2.ke);
T = double(vpasolve(3*x+ic1.m*besselk(1,ic1.m/x)/besselk(2,ic1.m/x)+...
    3*x+ic2.m*besselk(1,ic2.m/x)/besselk(2,ic2.m/x)==ic1.e+ic2.e,x,temp));
clear x
% MJ func
MJDF=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
    .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);
%% time
t=0; 
step_inc=1; step=step_inc;  idx=1;

%%  Other properties 

f_max=0.001; % maximal frequancy
cs_p=0; % cross section dependencies
coll_num=0; coll_num_Mm=0; %collision count

%% histogram and data saving 
edges=10.^(-16:0.1:16); % histogram edges 
data.f_sim_1=zeros(max_time,length(edges)-1);
data.f_sim_2=zeros(max_time,length(edges)-1);

%% Random diraction
dir_mat = randdir_matrix(N);
ind_dir = 1;
%% time  
step_inc=1; step=step_inc; idx=0; trys=0;
t=0;
temp_coll=0; temp_try=0;
%% loop 
while t<max_time
    trys=trys+1;
    id1=randi(2*N);
    id2=randi(2*N);
    while id2==id1
        id2=randi(2*N);
    end

    if id1>N
        M1=ic2.m;
    else
        M1=ic1.m;
    end

    if id2>N
        M2=ic2.m;
    else
        M2=ic1.m;
    end
    % Assign energy and momentum to the colliding particles
    e1=Enb(id1);
    e2=Enb(id2);

    p1=sqrt((Enb(id1))^2-M1^2)*dir_mat(ind_dir,:);
    ind_dir=mod(ind_dir+1,N);
    if ind_dir==0
        ind_dir=1;
    end

    p2=sqrt((Enb(id2))^2-M2^2)*dir_mat(ind_dir+1,:);
    ind_dir=mod(ind_dir+1,N);
    if ind_dir==0
        ind_dir=1;
    end
    
    f_now =fparam*collision_freq(p1,e1,p2,e2,cs_p); % frequancy calculation 
    f_max = max([f_now f_max]); % maximum frequency tracking

    dt=1/(N*f_max);
    t=t+dt;
    % Colision
    if f_now>f_max*rand
        coll_num=coll_num+1; % count the collision number performed 
        if ((M1==ic1.m && M2==ic2.m)||(M1==ic2.m && M2==ic1.m))
            coll_num_Mm=coll_num_Mm+1;
        end
        n0=dir_mat(ind_dir+1,:); % choose the direction in the CM frame after the collision
        ind_dir=mod(ind_dir+1,N); % advancing the Direction matrix index.. 
        if ind_dir==0
            ind_dir=1;
        end

        [Enb(id1), ~,Enb(id2), ~] = RelCol(e1, p1,e2, p2, n0);
    end

    if t>step
        step=step+step_inc;
        idx=idx+1;

        [data.bins,data.f_sim_1(idx,:)] = histlog(Enb(1:N)-ic1.m,edges);
        [data.bins,data.f_sim_2(idx,:)] = histlog(Enb((N+1):(2*N))-ic2.m,edges);
    
        data.coll_num(idx)=coll_num;
        data.coll_num_Mm(idx)=coll_num_Mm;
        data.time(idx)=t;

        disp({'time step', 'coll_num' , 'trys';...
            t coll_num trys-temp_try; ...
            'f_max' 'dt' 'd coll' ; ...
            f_max dt coll_num-temp_coll})
        temp_coll=coll_num;
        temp_try=trys;
    end
end
data.Teo1=MJDF(data.bins,ic1.m,T);
data.Teo2=MJDF(data.bins,ic2.m,T);

folder='saves';
mkdir(folder);
str=['Mass_ratio_' num2str(ic2.m/ic1.m)...
    '_gb_' num2str(ic1.gb,'%.0e') ...
    '_f_param_' num2str(fparam,'%.0e')...
    '_Nop_' num2str(N) ...
    '_max_time_' num2str(max_time) '.mat'];

save(fullfile(folder,str),'-mat')
end