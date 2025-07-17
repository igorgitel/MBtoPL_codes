clear 
close all 
clc 

%-----------------------------------------------
% Simulation of a single-species relativistic gas
% using func_One_type_prep.m
%-----------------------------------------------

% --- Physical parameters ---
m = 1;                  % Rest mass of the particles
gb = 10;               % Reduced momentum (gamma * beta), dimensionless

% Compute Lorentz gamma factor and corresponding energy
b = sqrt(gb^2 / (1 + gb^2));    % Velocity as a fraction of c
g = 1 / sqrt(1 - b^2);          % Lorentz factor
e = m * g;                      % Total energy per particle
ke = e - m;                     % Kinetic energy

% --- Simulation parameters ---
N = 1e6;                        % Number of particles
step_inc = 1e-3;                % Time interval between diagnostics (snapshots)
max_time = 2;                   % Total simulation time

% --- Run the simulation ---
% This function simulates elastic binary collisions and tracks
% the energy distribution over time.
% Outputs:
% - Enb: energy array of all particles
% - data: diagnostics structure (histograms, mean energy, etc.)
% - ic: initial condition struct (mass, energy, particle count, etc.)
[Enb, data, ic] = func_One_type_prep(m, gb, N, step_inc, max_time);
%% Plotting 

length(data)
ind=[1e-3 5e-3 1.5e-2 5e-2 1e-1 3e-1 1]*1e3;  


FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

for j=ind
    x=data(j).bins_Ek;
    y=data(j).f_sim_Ek;
    plot(x,y,'o-', LineWidth=1, MarkerSize=4)
    hold on 
end
plot(data(j).bins_Ek,data(j).Teo_Ek,'-', LineWidth=1.5)
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=22;
xlabel('$E_k=E-m$',Interpreter='latex')
ylabel('$f(E_k)$',Interpreter='latex')
xlim([6e-3 1e2])
ylim([1e-6 1e0])



text(1e1,5e-1,'t=0','FontSize',14);

t=text(1.5e-1,4e-5,{['$' num2str(data(ind(1)).time,"%.3f") ',\,\,' num2str(data(ind(1)).coll_num,"%.3f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#0072BD";

t=text(6e-2,8e-5,{['$' num2str(data(ind(2)).time,"%.3f") ',\,\,' num2str(data(ind(2)).coll_num,"%.2f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#D95319";

t=text(2.5e-2,1.5e-4,{['$' num2str(data(ind(3)).time,"%.2f") ',\,\,' num2str(data(ind(3)).coll_num,"%.2f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#EDB120";

t=text(2e-2,2.7e-4,{['$' num2str(data(ind(4)).time,"%.2f") ',\,\,' num2str(data(ind(4)).coll_num,"%.1f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#7E2F8E";

t=text(1e-2,5e-4,{['$' num2str(data(ind(5)).time,"%.1f") ',\,\,' num2str(data(ind(5)).coll_num,"%.1f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#77AC30";

t=text(0.8e-2,8.5e-4,{['$' num2str(data(ind(6)).time,"%.1f") ',\,\,' num2str(data(ind(6)).coll_num,"%.1f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#4DBEEE";

t=text(1.5e-2,8e-3,{['$t=' num2str(data(ind(7)).time,"%.0f") ',\,\,N_{cpp}=' num2str(data(ind(7)).coll_num,"%.0f") '$']},'FontSize',14);
t.FontSize=14;
t.Interpreter='latex';
t.Color="#A2142F";

   x = [0.25 0.15];
   y = [0.8 0.56];
g=annotation('textarrow',x,y,'String','Maxwell-Juttner DF');
g.FontSize=14;
title("b)")

exportgraphics(fig1,['MJ_Evol_Ek__initial_gb_' num2str(gb) '_UR_v2.png'],'Resolution',300)

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig2=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

for j=ind
    x=data(j).bins_mu;
    y=data(j).f_sim_mu;
    plot(x,y,'o-', LineWidth=1, MarkerSize=4)
    hold on 
end
plot(data(j).bins_mu,data(j).Teo_mu,'-', LineWidth=1.5)
ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=22;

xlabel('$\mu$',Interpreter='latex')
ylabel('$f(\mu)$',Interpreter='latex') 
xlim([1e-1 1e2])
ylim([1e-6 1e0])

text(1e1,5e-1,'t=0','FontSize',14);
text(1e0,1e-5,['t=' num2str(ind(1).*step_inc)],'FontSize',14);
text(5e-1,2e-5,[num2str(ind(2).*step_inc)],'FontSize',14);
text(3e-1,4e-5,[num2str(ind(3).*step_inc)],'FontSize',14);
text(3.5e-1,1.1e-4,[num2str(ind(4).*step_inc)],'FontSize',14);
text(2e-1,8e-5,[num2str(ind(5).*step_inc)],'FontSize',14);
text(1.5e-1,1.3e-4,[num2str(ind(6).*step_inc)],'FontSize',14);
text(2e-1,0.8e-3,[num2str(ind(7).*step_inc)],'FontSize',14);
x = [0.25 0.15];
y = [0.8 0.42];
g=annotation('textarrow',x,y,'String','Maxwell-Juttner DF');
g.FontSize=14;
title("a)")
exportgraphics(fig2,['MJ_Evol_mu__initial_gb_' num2str(gb) '_UR_v2.png'],'Resolution',300)
