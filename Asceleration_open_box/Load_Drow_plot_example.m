clear
dir Saves/
colors={'#0072BD','#D95319','#EDB120', '#7E2F8E', '#77AC30', '#4DBEEE', '#A2142F', ...
    '#D8B195','#F67280','#C06C84','#6C5B7B','#355C7D'};
color_ind=1;

MBDF=@(e,T) 2/(sqrt(pi*T^3)).*exp(-e/T).*sqrt(e);


FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

% linetype={'-','--','-.',':','-','--','-.',':'};
linetype={'-','-','-','-','-','-','-','-''-','-'};

% subplot(2,2,1)
f_param=1e0;
num_par_inj=1e4;
max_time=2e3;
cs_p=0;

str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')];
load(fullfile("Saves",[str '_mat_dir.mat']));

ind=[1 10 50 250 450 850 2e3];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins/1e8,data(ind(k)).f_sim);
%     q0=trapz(data(ind(k)).bins/1e8,data(ind(k)).f_sim/q);
    h=stairs(data(ind(k)).bins/1e8,data(ind(k)).f_sim/q);
    h.Color=colors{color_ind};
    color_ind=mod(color_ind+1,length(colors));
    if color_ind==0 
        color_ind=1;
    end

    h.LineStyle='-';
    hold on
    h.LineWidth=1;
end

f_param=1e0;
num_par_inj=1e4;
max_time=1;
step_inc=0.1;
cs_p=0;

str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')];
S=load(fullfile("Saves",[str '_mat_dir.mat']),'data');
data01=S.data;
% [data01]=Real_time_func_Evolution_single_injection_dir_matrix(f_param,num_par_inj,max_time,cs_p,step_inc);
for k=1:length(ind)
    q=trapz(data01(1).bins/1e8,data01(1).f_sim);
    h=stairs(data01(1).bins/1e8,data01(1).f_sim/q);
    %     h.LineStyle=linetype(k);
    %     h.LineStyle='-';
    hold on
    h.Color=colors{color_ind};
    color_ind=mod(color_ind+1,length(colors));
    if color_ind==0 
        color_ind=1;
    end
    h.LineWidth=1;
%     xd = data(ind(k)).bins(data(ind(k)).f_sim~=0);
%     yd = data(ind(k)).f_sim(data(ind(k)).f_sim~=0)/q;
%     q=trapz(xd,xd.*yd);
%     i1=round(2/3*length(xd));
%     text(xd(i1),yd(i1)*5,num2str(data(ind(k)).time),'FontSize',16,'fontweight','bold',Interpreter='latex');
%     data(ind(k)).time
end
% plot(data(ind(k)).bins,MBDF(data(ind(k)).bins,T),'-')
xlabel('$E_k/mc^2$','Interpreter','latex')
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
text('FontSize',16,'Interpreter','latex','String','$t = 0.1$',...
    'Position',[0.777855413904808/1e8 7.04682699434797*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','1',...
    'Position',[2.96596799051017/1e8 0.238176590947677*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','10',...
    'Position',[68.0695881005449/1e8 0.0136201602456586*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','50',...
    'Position',[1023.08700581328/1e8 0.000632574037147515*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','250',...
    'Position',[16198.5337363002/1e8 3.47435681857438e-05*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','450',...
    'Position',[220.962192577822/1e8 2.03006156512141e-06*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','850',...
    'Position',[461.701815350814/1e8 4.76929316072696e-07*Norm_temp 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String',{'$t_{eq}=2000$'},...
    'Position',[731.379756357976/1e8 9.02483429972361e-08*Norm_temp 0]);


exportgraphics(fig1,'Clasic_fermi_asc_1.jpg','Resolution',300)
%---------------
% subplot(2,2,2)

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig2=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

f_param=1e0;
num_par_inj=1e5;
max_time=1e1;
step_inc=0.01;
cs_p=0;


str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')...
    '_step_inc_' num2str(step_inc,'%.2f')];
load(fullfile('Saves',[str '_mat_dir.mat']))


ind=[1e1 1e2 3e2 1e3];

x=[1 1 1 1.5e-1]; % the shift
y=[2 5e-1 8e-1 6];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    h=stairs(x(k)*data(ind(k)).bins/1e8,y(k)*data(ind(k)).f_sim/(q));
    h.LineStyle=linetype(k);
    hold on
    h.LineWidth=1;
%     xd = data(ind(k)).bins(data(ind(k)).f_sim~=0);
%     yd = data(ind(k)).f_sim(data(ind(k)).f_sim~=0)/q;
%     [maxi, i1]=max(yd);
%     text(xd(i1)*2,yd(i1)*2,num2str(data(ind(k)).time),'FontSize',16,'fontweight','bold',,Interpreter='latex');
    data(ind(k)).time
end

ylim([1e-6 3e1])
xlim([6e-11 5e-6])
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
text('FontSize',16,'Interpreter','latex','String','$t=0.1$',...
    'Position',[3.95651323484695/1e8 2.41803034897478e-05 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','$1$',...
    'Position',[38/1e8 0.000277003347634997 0]);

% Create arrow
annotation(fig2,'textarrow',[0.65 0.62],[0.69 0.62],...
    'FontSize',16,'Interpreter','latex','String','$3$');

% Create arrow
annotation(fig2,'textarrow',[0.76 0.72],[0.61 0.54],...
    'FontSize',16,'Interpreter','latex','String','$10$');
title('a)')
exportgraphics(fig2,'Clasic_fermi_asc_2.jpg','Resolution',300)

%----------------------------------------
% subplot(2,2,3)
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig3=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);
f_param=1e0;
num_par_inj=1e4;
max_time=2e3;
cs_p=0;

str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')];
load(fullfile("Saves",[str '_mat_dir.mat']));

ind=[10 50 250 450];

x=[1300 60 3 1]; % the shift
% y=[5 1 1 10];
clear leg
for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    % h=stairs(x(k)*data(ind(k)).bins/1e8,data(ind(k)).f_sim/(q*x(k)));
    h=plot(x(k)*data(ind(k)).bins/1e8,data(ind(k)).f_sim/(q*x(k)));
    h.LineStyle=linetype(k);
    hold on
    h.LineWidth=1;
    data(ind(k)).time
    leg{k}=num2str(data(ind(k)).time);
end

xlim([5e-6 3e-2])
ylim([8e-10 5e-5])
xticks([10.^(-16:1:16)])
yticks([10.^(-16:1:16)])
xlabel('$x$','Interpreter','latex')
ylabel('$f(x)$',Interpreter='latex')
plg=legend(leg,'Location','southwest',Interpreter="latex");
title(plg,'t');
ax=gca;
ax.XScale='log';
ax.YScale='log';
t=text(1e5/Norm_temp,1e-5,'$\propto \sqrt{x}e^{-\sqrt{x}}$','Interpreter','latex','FontSize',16,'fontweight','bold');
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

ind=[450 853 2000];

x=[3 1.5 1];

for k=1:length(ind)
    q=trapz(data(ind(k)).bins,data(ind(k)).f_sim);
    h=stairs(x(k).*data(ind(k)).bins/1e8,data(ind(k)).f_sim/(x(k)*q));
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
    'Position',[1000/1e8 8e-07 0]);

% Create text
text('FontSize',16,'Interpreter','latex','String','$t_{eq}=2000$',...
    'Position',[885/1e8 4.2e-07 0]);
% Create arrow
annotation(fig4,'textarrow',[0.25 0.20],[0.45 0.58],...
    'FontSize',16,'Interpreter','latex','String','MB');
title('c)')
exportgraphics(fig4,'Clasic_fermi_asc_4.jpg','Resolution',300)
%% Mean ------------------


f_param=1e0;
num_par_inj=1e4;
max_time=2e3;
cs_p=0;

str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')];
load(fullfile("Saves",[str '_mat_dir.mat']));

for k=1:length(data)
    time(k)=data(k).time;
    kin_energy(k)=data(k).kin_energy;
end

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig5=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

h=stairs(time,kin_energy/1e8,'-');
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

ylabel('$\langle E_k \rangle/mc^2$',Interpreter='latex')
xlabel('$t$',Interpreter='latex')

xticks(10.^(0:4))
yticks(10.^(-8:4))

annotation(fig5,'textarrow',[0.65 0.78],[0.80 0.82],...
    'FontSize',16,'Interpreter','latex','String','$\langle E_{k} \rangle \propto t^2 $');
ax.FontSize=18;
title(['d)'])
exportgraphics(fig5,'Clasic_fermi_asc_5.jpg','Resolution',300)
%% Cumulative DF Constant injection DF 

f_param=1e0;
num_par_inj=1e4;
max_time=2e3;
cs_p=0;

str=['Real_time_Single_injection_number_of_particles_injected_' num2str(num_par_inj,'%1.1e')...
    '_cs_power_' num2str(cs_p,'%1.0f') ...
    '_maximum_time_steps_' num2str(max_time,'%1.0e') ...
    '_relative_frequancy_parameter_' num2str(f_param,'%1.0e')];
load(fullfile("Saves",[str '_mat_dir.mat']));
data(1);

bins=data(1).bins;
DFCI=data(1).f_sim;

for i=2:length(data)
    DFCI=DFCI+data(i).f_sim;
end

for j=1:1e4
    indx=randi([1995 2000]);
    DFCI=DFCI+data(indx).f_sim;
end
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig6=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

sb1=subplot(4,1,[1 3]);
sb1.Position=[0.1300    0.351    0.7750    0.5959];
q=trapz(bins,DFCI)
h=stairs(bins/1e8,DFCI/q);
h.LineWidth=1;
h.LineStyle='-';

ax=gca;
ax.XScale='log';
ax.YScale='log';
xl=[1e-1 1e7]./1e8;
xlim(xl)
yl=[1e1 1e5]./1e8;
ylim(yl)

xticks(10.^[(log10(xl(1))):(log10(xl(2)))])
xticklabels([])
yticks(10.^[(log10(yl(1))+1):1:(log10(xl(2))-1)])

xlabel('$E_k$',Interpreter='latex')
ylabel('$f(E_k)$',Interpreter='latex')
ax.FontSize=18;
grid off

x=bins; y=DFCI;
sb2=subplot(4,1,4);
sb2.Position=[0.1300    0.1600    0.7750    0.19];
% sb2.OuterPosition=[0    0.02   1.0000    0.335];
h=3;
for i=(h+1):(length(y))
    dy(i)=log10(y(i))-log10(y(i-h));
    dy(i)=dy(i)/(log10(x(i))-log10(x(i-h)));
end
h=stairs(bins/1e8,dy,'LineWidth',1.5);
h.MarkerSize=3;
ax=gca;
ax.XScale='log';
ylim([-1 0.2])
xlim(xl)
yline(-0.5,LineWidth=1)
xlabel('$E_k/mc^2$',Interpreter='latex')
ylabel('PL index',Interpreter='latex')
ax.FontSize=18;
xticks(10.^[(log10(xl(1))):(log10(xl(2)))])
yticks([-1 -0.5 0])

exportgraphics(fig6,'Clasic_fermi_asc_6.jpg','Resolution',300)