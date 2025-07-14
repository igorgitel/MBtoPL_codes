function [data ,str,ic,icb] = SPI_RED_cascade(ic,icb,...
    maxnumofpar,numofpar,...
    max_time,folder)
addpath('Basic functions/')
%% Initial conditions
icb.g=sqrt(icb.gb^2+1);
icb.e=icb.m*icb.g;
icb.ek=icb.e-icb.m;

% Energy array
ic.b=sqrt(ic.gb^2/(1+ic.gb^2));
ic.g=sqrt(ic.gb^2+1);
ic.e=ic.m*ic.g;
ic.ek=ic.e-ic.m;

En=zeros(1,10*maxnumofpar); % array def
En(1:numofpar)=ic.e*ones(1,numofpar); % The injection

% redusction options
num_of_par_after_red=maxnumofpar/10; %number of particle after reduction;
numofpar_initially=numofpar; % for saving

DF_factor=1; % factor for the distribtution fanction after reduction

% background
[Enb,~,~]=func_One_type_prep(icb.m,icb.gb,icb.N,icb.folder);
drawnow

% time
t=0; step_inc=0.01; step=step_inc; idx=0; %max_time=1e4;
N=numofpar;

% frequncies
f_max=0;
f_param=1e0;

% Intermidiate output
temp1=0; cs_p=0;

% min energy stop option
ratio_nop_below_min_energy=0;

% edges
edges=10.^(-16:0.1:16);


% The structure
data.bins=zeros(1,length(edges)-1);
data.fsim=zeros(max_time/step_inc,length(edges)-1);
data.fsim10=zeros(max_time/step_inc,length(edges)-1);
data.fsim50=zeros(max_time/step_inc,length(edges)-1);
data.fsim100=zeros(max_time/step_inc,length(edges)-1);
data.EffParNum=zeros(1,max_time/step_inc);
data.time=zeros(1,max_time/step_inc);
data.DF_factor=zeros(1,max_time/step_inc);
data.kin_energy=zeros(1,max_time/step_inc);
erff=zeros(1,1e6);

% random diractions Matrix
Max_dir_ind=1e6;
dir_mat = randdir_matrix(Max_dir_ind);
ind_dir=1;

%%
numberofcollision=0;

while t<max_time
    id1 = randi(numofpar);
    id2 = randi(icb.N);

    % colision prepare
    e1=En(id1);
    e2=Enb(id2);
    p1=sqrt((En(id1))^2-ic.m^2)*dir_mat(ind_dir,:);
    ind_dir=ind_dir+1;
    if ind_dir>=Max_dir_ind
        ind_dir=1;
    end
    p2=sqrt((Enb(id2))^2-icb.m^2)*dir_mat(ind_dir+1,:);
    ind_dir=ind_dir+1;
    if ind_dir>=Max_dir_ind
        ind_dir=1;
    end

    f_now=collision_freq(p1,e1,p2,e2,f_param,cs_p);
    f_max=max([f_now f_max]);
    dt=(1/N)*(1/f_max);
    t=t+dt;

    if f_now>f_max*rand
        numberofcollision=numberofcollision+1; %number of colision
        te=e1+e2;

        n0=dir_mat(ind_dir+1,:);
        ind_dir=ind_dir+1;
        if ind_dir>=Max_dir_ind
            ind_dir=1;
        end

        [ea1, ~,ea2, ~] = RelCol(e1, p1,e2, p2, n0);

        erff(numberofcollision)=(te-(ea1+ea2))/e2;

        numofpar=numofpar+1;
        En(id1)=ea1; En(numofpar)=ea2;
    end

    if t>step
        idx=idx+1;

        kin_energ=mean(En(1:numofpar)-ic.m);

        % minimal energy criterial
        kin_energy_array=En(1:numofpar)-ic.m;
        nop_below_min_energy=sum(kin_energy_array<10*icb.ke);
        ratio_nop_below_min_energy=nop_below_min_energy/numofpar;


        % histogram and data
        [data.bins,data.fsim10(idx,:),data.fsim50(idx,:),...
            data.fsim100(idx,:),data.fsim(idx,:),~] = ...
            low_err_histlog(En-ic.m,edges);
% 

        % DF factor
        data.fsim(idx,:)=DF_factor*data.fsim(idx,:);
        data.fsim10(idx,:)=DF_factor*data.fsim10(idx,:);
        data.fsim50(idx,:)=DF_factor*data.fsim50(idx,:);
        data.fsim100(idx,:)=DF_factor*data.fsim100(idx,:);

        data.DF_factor(idx)=DF_factor;
        data.EffParNum(idx)=DF_factor.*numofpar;
        data.kin_energy(idx)=kin_energ;
        data.time(idx)=t;

        N=numofpar;

        % reduction mechnism
        if numofpar>maxnumofpar
            red_factor=numofpar/num_of_par_after_red;
            ind=randi([1 numofpar],1,num_of_par_after_red);

            En(1:num_of_par_after_red)=En(ind);
            En((num_of_par_after_red+1):end)=0;

            numofpar=num_of_par_after_red;
            DF_factor=DF_factor*red_factor;
            N=numofpar;

        end


        disp({'time step' 'inc particle' 'kin energy' 'dt'; ...
            step DF_factor.*(numofpar-temp1)  kin_energ dt;...
            'B_kin_energy' 'Eff_num_par'  'num. col.' 'f_max'; ...
            icb.ke data.EffParNum(idx) numberofcollision f_max ; ...
            'nop_bme'  'r_nop_bme' 'DF_factor' 'gg'; ...
            nop_below_min_energy ratio_nop_below_min_energy DF_factor 'g'})

        temp1=numofpar;
        step=step+step_inc;
    end
end
%%
% figure
% histogram(erff,10.^(-16:0.1:1));
% ax=gca;
% ax.XScale='log';
% ax.YScale='log';
% title(num2str(mean(erff),'%1.1e'))
%%
str=['init_gb_' num2str(ic.gb,'%1.0e') ...
    '_par_num_' num2str(numofpar_initially,'%1.0e') ...
    '_max_num_' num2str(maxnumofpar,'%1.0e') ...
    '_cs_p_' num2str(cs_p) ...
    '_red_rand_dir'];
% folder='RD_saves_05';
mkdir(folder)
save(fullfile(folder,[str '.mat']),"data","ic","str",'icb','-v7.3')

toc
%  evol_Mat_anim(data,str,folder)
%  anima_s(s,str,ic)
end