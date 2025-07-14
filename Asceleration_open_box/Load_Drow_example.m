clear
% addpath('Basic functions\')
folder='saves_v6';
dan=dir(folder)
colors={'#0072BD','#D95319','#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', ...
    '#D8B195','#F67280','#C06C84','#6C5B7B','#355C7D'};
color_ind=1;

MBDF=@(e,T) 2/(sqrt(pi*T^3)).*exp(-e/T).*sqrt(e);



% linetype={'-','--','-.',':','-','--','-.',':'};
linetype={'-','-','-','-','-','-','-','-''-','-'};

% subplot(2,2,1)
% f_param=1e0;
% num_par_inj=1e5;
% max_time=2.5e3;
% step_inc=0.05;
% cs_p=0;
% 
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% 
% load(fullfile(folder,[str '_mat_dir.mat']));
in=5;
load(fullfile(dan(in).folder,dan(in).name))
%%
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);
x=data(1).bins;
y=data(1).f_sim;
y=y./trapz(x,y);
e_first=trapz(x,x.*y);
icb.m=1e6;
icb.e=icb.m/2;
ic.e=0.5;
Tb=2/3*icb.e;
Ta=2/3*ic.e;


ind=[1 2 [1 10 50 250 450 850]*20+1 2.5e3*20];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins/1e8,data(ind(k)).f_sim);
%     q0=trapz(data(ind(k)).bins/1e8,data(ind(k)).f_sim/q);
    h=plot(data(ind(k)).bins/1e8,data(ind(k)).f_sim/q);
    h.Color=colors{color_ind};
    color_ind=mod(color_ind+1,length(colors));
    if color_ind==0 
        color_ind=1;
    end

    h.LineStyle='-';
    hold on
    h.LineWidth=1;
end
x=data(ind(k)).bins;

q=1e-8;
h=plot(x/1e8,MBDF(x,Tb)/q,'--',Color='k');
q=1e-8;
h=plot(x/1e8,MBDF(x,Ta)/q,'--',Color='k');

% f_param=1e0;
% num_par_inj=1e4;
% max_time=1;
% step_inc=0.1;
% cs_p=0;
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% 
% danie=load(fullfile(folder,[str '_mat_dir.mat']));
% data01=danie.data;
% clear danie
% 
% [data01]=Real_time_func_Evolution_single_injection_dir_matrix(f_param,num_par_inj,max_time,cs_p,step_inc);

% for k=1:length(ind)
%     q=trapz(data01(1).bins/1e8,data01(1).f_sim);
%     h=plot(data01(1).bins/1e8,data01(1).f_sim/q);
%     %     h.LineStyle=linetype(k);
%     %     h.LineStyle='-';
%     hold on
%     h.Color=colors{color_ind};
%     color_ind=mod(color_ind+1,length(colors));
%     if color_ind==0 
%         color_ind=1;
%     end
%     h.LineWidth=1;
%     xd = data(ind(k)).bins(data(ind(k)).f_sim~=0);
%     yd = data(ind(k)).f_sim(data(ind(k)).f_sim~=0)/q;
%     q=trapz(xd,xd.*yd);
%     i1=round(2/3*length(xd));
%     text(xd(i1),yd(i1)*5,num2str(data(ind(k)).time),'FontSize',16,'fontweight','bold',Interpreter='latex');
%     data(ind(k)).time
% end
% plot(data(ind(k)).bins,MBDF(data(ind(k)).bins,T),'-')
xlabel('$E_k$','Interpreter','latex')
ylabel('$f(E_k)$',Interpreter='latex')
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=18;
grid off
% 

xticks(10.^(-10:2:1))
yticks(10.^(-9:2:9))

T=(1e6+1)/3;
Norm_temp=1e8;
ylim([1e-9 5e1]*Norm_temp)
xlim([1e-10 5e-2])

% Create text
text('FontSize',16,'Interpreter','latex','String',['$t=' num2str(data(ind(1)).time,'%0.0f') '$'],...
    'Position',[6e-9 1e5 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$t=' num2str(data(ind(2)).time,'%0.2f') '$'],...
    'Position',[0.777855413904808/1e8 3.04682699434797*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(3)).time,'%0.0f') '$'],...
    'Position',[2.96596799051017/1e8 0.238176590947677*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(4)).time,'%0.0f') '$'],...
    'Position',[68.0695881005449/1e8 0.0136201602456586*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(5)).time,'%0.0f') '$'],...
    'Position',[1023.08700581328/1e8 0.000632574037147515*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(6)).time,'%0.0f') '$'],...
    'Position',[16198.5337363002/1e8 3.47435681857438e-05*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(7)).time,'%0.0f') '$'],...
    'Position',[120.962192577822/1e8 2.03006156512141e-06*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(8)).time,'%0.0f') '$'],...
    'Position',[251.701815350814/1e8 4.76929316072696e-07*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$t_{eq}=' num2str(data(ind(9)).time,'%0.0f') '$'],...
    'Position',[1531.379756357976/1e8 9.02483429972361e-08*Norm_temp 0]);


exportgraphics(fig1,'Clasic_fermi_asc_1.jpg','Resolution',300)
%%

%---------------
% subplot(2,2,2)

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig2=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

% f_param=1e0;
% num_par_inj=1e4;
% max_time=1e1;
% step_inc=0.01;
% cs_p=0;
% 
% 
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% % [~]=Real_time_func_Evolution_single_injection_dir_matrix(f_param,num_par_inj,max_time,cs_p,step_inc);
% 
% load(fullfile(folder,[str '_mat_dir.mat']))

% ic.m=1;
% ic.e=ic.m/2;
% ic.N=1e5;
% 
% str=['mass_' num2str(ic.m,'%1.0e')...
%      '_energy_' num2str(ic.e,'%1.0e')...
%      '_particle_number_' num2str(ic.N,'%1.0e')];
% 
% Array=load(fullfile('saves_one_type',[str '.mat']));

ind=[1 5e0 2e1+1 5e1+1 2e2+1];

% x=[1 1 1 1.5e-1]; % the shift for initial delta
% y=[2 5e-1 8e-1 6];

% loglog(Array.data(end).bins/1e8, Array.data(end).f_sim,'-',LineWidth=1)
% bins=data(1).bins;
% loglog(bins/1e8,MBDF(bins,Ta))
% hold on

x=[1 1 0.5 0.3 0.04]; % the shift for initial MB 
y=[1 1 2 5 100];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    h=plot(x(k)*data(ind(k)).bins/1e8,y(k)*data(ind(k)).f_sim/(q));
    h.LineStyle=linetype(k);
    hold on
    h.LineWidth=1;
%     xd = data(ind(k)).bins(data(ind(k)).f_sim~=0);
%     yd = data(ind(k)).f_sim(data(ind(k)).f_sim~=0)/q;
%     [maxi, i1]=max(yd);
%     text(xd(i1)*2,yd(i1)*2,num2str(data(ind(k)).time),'FontSize',16,'fontweight','bold',,Interpreter='latex');
    data(ind(k)).time
end

ylim([1e-7 3e1])
xlim([5e-11 5e-6])
xticks(10.^(-16:2:16))
yticks(10.^(-16:2:16))
xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$',Interpreter='latex')
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=18;
grid off

% Create text
text('FontSize',16,'Interpreter','latex','String',['$t=' num2str(data(ind(1)).time,'%0.0f') '$'],...
    'Position',[2e-8 1e-4 0]);
% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(2)).time,'%0.1f') '$'],...
    'Position',[9e-8 3e-5 0]);
% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(3)).time,'%0.0f') '$'],...
    'Position',[3e-7 2e-6 0]);
% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(4)).time,'%0.1f') '$'],...
    'Position',[5e-7 1e-6 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',['$' num2str(data(ind(5)).time,'%0.1f') '$'],...
    'Position',[4e-7 1e-3 0]);

title('a)')
exportgraphics(fig2,'Clasic_fermi_asc_2.jpg','Resolution',300)

%%

%----------------------------------------
% subplot(2,2,3)
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig3=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

% f_param=1e0;
% num_par_inj=1e4;
% max_time=2.3e3;
% step_inc=1;
% cs_p=0;
% 
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% load(fullfile(folder,[str '_mat_dir.mat']));

ind=[10 50 250 450]*20+1;

x=[1500 60 3 1]; % the shift
% y=[5 1 1 10];
clear leg
for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    h=plot(x(k)*data(ind(k)).bins/1e8,data(ind(k)).f_sim/(q*x(k)));
    % cftool(log10(data(850).bins),log10(data(850).f_sim))
    % log10(A)+a*x-log10(exp(1))*10^(b*(x-log(T)))
    h.LineStyle=linetype(k);
    hold on
    h.LineWidth=1;
    data(ind(k)).time
    leg{k}=num2str(data(ind(k)).time);
end

xlim([1e-6 3e-2])
ylim([1e-10 5e-5])
xticks([10.^(-16:1:16)])
yticks([10.^(-16:1:16)])
xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$',Interpreter='latex')
plg=legend(leg,'Location','southwest',Interpreter="latex");
title(plg,'t');
ax=gca;
ax.XScale='log';
ax.YScale='log';
t=text(1e5/Norm_temp,1e-5,'$\propto \sqrt{x}e^{\sqrt{x}}$','Interpreter','latex','FontSize',16,'fontweight','bold');
ax.FontSize=18;
grid off
title('b)')
exportgraphics(fig3,'Clasic_fermi_asc_3.jpg','Resolution',300)
%----------------------------------------
% subplot(2,2,4)
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig4=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

ind=[[450 853]*20+1 2.5e3*20];

x=[3 1.5 1];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    h=plot(x(k).*data(ind(k)).bins/1e8,data(ind(k)).f_sim/(x(k)*q));
    h.LineStyle=linetype(k);
    hold on
    h.LineWidth=1;
%     xd = data(ind(k)).bins(data(ind(k)).f_sim~=0);
%     yd = data(ind(k)).f_sim(data(ind(k)).f_sim~=0)/q;
%     [maxi, i1]=max(yd);
%     t=text(xd(i1)*2,yd(i1)*2,num2str(data(ind(k)).time),'FontSize',16,'fontweight',,);

    data(ind(k)).time
end
% plot(data(ind(k)).bins,MBDF(data(ind(k)).bins,T),'-')
xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$',Interpreter='latex')
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=18;
grid off
% 
% title('The distributions in energy')


xticks([10.^(-4:1:9)]./1e8)
yticks(10.^(-10:-4))
T=(1e6+1)/3;
% T=T-2.8e5
q=trapz(data(ind(k)).bins,MBDF(data(ind(k)).bins,T));
h=plot(data(ind(k)).bins/1e8,MBDF(data(ind(k)).bins,T),':');


h.LineWidth=1;
ylim([1e-9 1e-5])
xlim([6e2 5e6]./1e8)

% Create text
text('FontSize',16,'Interpreter','latex','String','$t=450$',...
    'Position',[803/1e8 2.5e-06 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','850',...
    'Position',[2e-5 5e-07 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','$t_{eq}=2500$',...
    'Position',[5e-05 3e-07 0]);
% Create arrow
annotation(fig4,'textarrow',[0.20 0.15],[0.40 0.57],...
    'FontSize',16,'Interpreter','latex','String','MB');
title('c)')
exportgraphics(fig4,'Clasic_fermi_asc_4.jpg','Resolution',300)
%% Mean ------------------


% f_param=1e0;
% num_par_inj=1e4;
% max_time=2.3e3;
% step_inc=1;
% cs_p=0;
% 
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% load(fullfile(folder,[str '_mat_dir.mat']));

for k=1:length(data)
    time(k)=data(k).time;
    kin_energy(k)=data(k).kin_energy;
end

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig5=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

h=plot(time,kin_energy/1e8,'-');
h.LineWidth=1;

[xData, yData] = prepareCurveData( log10(time),log10(kin_energy/1e8) );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = (xData < log10(10)) | (xData > log10(200));
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
hold on
x=10.^(0:0.1:4);
h=plot(x,10.^(fitresult.p2)*(x).^(fitresult.p1),'--');
h.LineWidth=1;

ax=gca;
ax.XScale='log';
ax.YScale='log';
grid off


xlim([1e0 2e3])

ylabel('$\langle E_k \rangle$',Interpreter='latex')
xlabel('$t$',Interpreter='latex')

xticks(10.^(0:4))
yticks(10.^(-8:4))

annotation(fig5,'textarrow',[0.65 0.78],[0.80 0.82],...
    'FontSize',16,'Interpreter','latex','String','$\langle E_{k} \rangle \propto t^2 $');
ax.FontSize=18;
title(['d)'])
exportgraphics(fig5,'Clasic_fermi_asc_5.jpg','Resolution',300)
%% Cumulative DF Constant injection DF 

% f_param=1e0;
% num_par_inj=1e4;
% max_time=2.3e3;
% step_inc=1;
% cs_p=0;
% 
% str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
%     '_cs_power_' num2str(cs_p,'%1.0f') ...
%     '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
%     '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
%     '_step_inc_' num2str(step_inc,'%.2f')];
% load(fullfile(folder,[str '_mat_dir.mat']));

MBDF=@(e,T) 2/(sqrt(pi*T^3)).*exp(-e/T).*sqrt(e);

Tb=2/3*icb.e
Ta=2/3*ic.e
bins=data(1).bins;

DFCI=0*MBDF(bins,Ta);
for i=1:length(data)
    q=trapz(data(i).bins,data(i).f_sim);
    DFCI=DFCI+data(i).f_sim/q;
end
DFCI=DFCI+1e5*MBDF(bins,Tb);

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig6=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

sb1=subplot(4,1,[1 3]);
sb1.Position=[0.1300    0.42    0.7750    0.58];
q=trapz(bins,DFCI);
h=plot(bins/1e8,DFCI/q);
hold on
h.LineWidth=1;
h.LineStyle='-';
q=5e3;
h=plot(bins/1e8,MBDF(bins,Ta)/q);
h.LineWidth=1;
h.LineStyle='--';

q=8.5e-1;
h=plot(bins/1e8,MBDF(bins,Tb)*q);
h.LineWidth=1;
h.LineStyle='--';


text(2e-9,1e-6,'Initial MB', Interpreter='latex')

text(2e-6,5e-7,'Background MB',Interpreter='latex')

ax=gca;
ax.XScale='log';
ax.YScale='log';
xl=[1e-3 1e7]./1e8;
xlim(xl)
yl=[1e-7 1e-3];
ylim(yl)

xticks(10.^[(log10(xl(1))):1:(log10(xl(2)))])
xticklabels([])
yticks(10.^[(log10(yl(1))+1):1:(log10(yl(2))-1)])

xlabel('$E_k$',Interpreter='latex')
ylabel('$F(E_k)$',Interpreter='latex')
ax.FontSize=18;
grid off

x=bins; y=DFCI;
sb2=subplot(4,1,4);
sb2.Position=[0.1300    0.22    0.7750    0.2];
% sb2.OuterPosition=[0    0.02   1.0000    0.335];
h=3;
for i=(h+1):(length(y))
    dy(i)=log10(y(i))-log10(y(i-h));
    dy(i)=dy(i)/(log10(x(i))-log10(x(i-h)));
end
h=plot(bins/1e8,dy);
h.LineWidth=1;
h.Marker='.';
h.MarkerSize=5;
ax=gca;
ax.XScale='log';
ylim([-1 0.2])
xlim(xl)
yline(-0.5,LineWidth=1)
xlabel(['$E_k$'],Interpreter='latex')
ylabel('PL index',Interpreter='latex')
ax.FontSize=18;
xticks(10.^[(log10(xl(1))):2:(log10(xl(2)))])
yticks([-1 -0.5 0])

exportgraphics(fig6,'Clasic_fermi_asc_6.jpg','Resolution',300)