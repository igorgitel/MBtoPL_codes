% Define simulation metadata to locate the file
N       = 1e5;       % Number of Type 1 particles
beta2   = 0.0995;       % Velocity (v/c) of background scatterers
E0      = 100;       % Initial energy of particles (e.g., gamma * m = 1.0)
stepinc = 1.0;       % Time between histogram snapshots

% Construct filename string (must match the one used in save block)
folder   = 'saves_two_type';  % Folder where results are stored
filename = ['TwoType_N_' num2str(N, '%.0e') ...
            '_beta2_' num2str(beta2, '%.0e') ...
            '_E0_' num2str(E0, '%.0e') ...
            '_step_' num2str(stepinc, '%.1f') '.mat'];

% Full path to file
filepath = fullfile(folder, filename);

% Load output data and configuration used in the simulation
load(filepath, 's', 'ic', 'edges', 'time_s');

% --------------------------------------------------------
% s        — Struct array with results per time step
% ic       — Input parameters: .p01 (initial momentum), .beta2 (target velocity)
% edges    — Histogram binning for energy distribution
% time_s   — Time control struct (e.g., stepinc)
% --------------------------------------------------------

for i=1:length(s)
%   meanEnerg1(i)=trapz(s(i).EnergyBE1,s(i).EnergyBE1.*s(i).DFE1);
    meanEnerg(i)=s(i).Energ_mean;
    time(i)=s(i).time;
    colnum(i)=s(i).coldone;
end

% Default color order in MATLAB (since R2014b)
defaultColors = [
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.4940, 0.1840, 0.5560;  % purple
    0.4660, 0.6740, 0.1880;  % green
    0.3010, 0.7450, 0.9330;  % light blue (cyan)
    0.6350, 0.0780, 0.1840   % dark red
];

color_index=1;
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

bins=s(1).EnergyBE1;
dfs=s(1).DFE1;

for i=2:length(s)
    dfs=dfs+s(i).DFE1;
end

ind=(6:2:11) * 1e3

sb1=subplot(4,1,[1 3]);
sb1.Position=[0.1300    0.33    0.7750    0.5959];
q=trapz(bins,dfs)
h=plot(bins,dfs,'Color', defaultColors(color_index,:),'DisplayName','Steady-state');
color_index=color_index+1;

h.LineWidth=1;
h.Marker='o';
h.MarkerSize=3;

ax=gca;
ax.XScale='log';
ax.YScale='log';
ax.FontSize=18;


% xticks(10.^(-1:1:3))
% yticks(10.^(0:1:6))

xlabel('$E_k$',Interpreter='latex')
ylabel('$F(E_k)$',Interpreter='latex')

xl=[1e0 1e20];
xlim(xl)
yl=[1e-20 1e2];
ylim(yl)

title(ax,'a)')

xticks(10.^[(log10(xl(1))):4:(log10(xl(2)))])
xticklabels([])
yticks(10.^[(log10(yl(1))+1):5:(log10(yl(2)))])
for in=ind
    x = s(in).EnergyBE1;
    f_x = s(in).DFE1; 

    hold on
    sh_str = sprintf('%.2e', double(s(in).Energ_mean));  % display only, keep double for label
    plot(x, f_x, '.','Color', defaultColors(color_index,:),'DisplayName', ['$E_0 = ' sh_str '$']);
    color_index = color_index+1;
    
end

legend(Location='southwest',Interpreter='latex', FontSize=14)
% subplot(4,2end,8)
sb2=subplot(4,1,4);
sb2.Position=[0.1300    0.1500    0.7750    0.18];


x=bins; y=dfs;
h=4;
for i=(h+1):(length(y))
    dy(i)=log10(y(i))-log10(y(i-h));
    dy(i)=dy(i)/(log10(x(i))-log10(x(i-h)));
end
h=plot(x,dy,'o-','LineWidth',2);
h.MarkerSize=3;
yline(-1)
xlabel('$E_k$', 'Interpreter', 'latex', 'FontSize',12);
ylabel('PL index')

% Limits
ylim([-1.5 0])
yticks([-1.5 -1.0 -0.5])

xlim(xl)
xticks(10.^[(log10(xl(1))):4:(log10(xl(2)))])

ax=gca;
ax.XScale='log';
ax.FontSize=18;

exportgraphics(fig1,'Ultra_reltivistic_evolution.jpg','Resolution',300)

%% 
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

% Axes styling
ax = gca;
% ax.FontName = 'Times New Roman';         % APS uses Times
ax.FontSize = 18;                          % APS recommends 8 pt font
% ax.LineWidth = 0.8;                       % Thin but visible
% ax.TickDir = 'out';
% ax.Box = 'off';                           % No top/right frame
% ax.XMinorTick = 'on';
% ax.YMinorTick = 'on';
ax.Box = 'on';
ax.TickDir = 'in';  % Optional but cleaner

hold on
legendEntries = {};
ind=(6:2:11) * 1e3
sh=ind*0+1;
k=1;
for i = ind
   E0(k)=s(i).Energ_mean;
   sh(k)=E0(k);
   k=k+1;
end

% sh(1)=2.48e26;
% sh(2)=1.11e29;
% sh(3)=1.24e34;

% shiftsx = 10.^([-5 0 0]);
% shiftsy = 10.^([3 0 0]);

k=1;
    color_index = 2;
for i = ind
    x = s(i).EnergyBE1;
    f_x = s(i).DFE1;

    % s_x=shiftsx(k);
    % s_y=shiftsy(k);

    x_scaled = x/sh(k);
    f_x_scaled = f_x .* sh(k) .* exp(-x_scaled);
    % for j = 1:length(f_x)        
    %     if f_x(j) == 0
    %         f_x_scaled(j)=0;
    %     else
    %         f_x_scaled(j) = log(f_x(j)) +log(sh(k)) - x(j) /sh(k);
    %     end
    % end
    mask=f_x_scaled ~= 0;

    
    % prod(f_x_scaled>=0)
    % double(min(f_x_scaled))
    title(ax,'b)')
    plot(x_scaled(mask), f_x_scaled(mask), '.','Color', defaultColors(color_index,:));
    color_index = color_index+1;
    x1=double(log10(x_scaled(mask)));
    y1=double(log10(f_x_scaled(mask)));

    plotFlag = false;
    [params, y_fit] = fitCustomModel(x1, y1, plotFlag);
    %params
    % cftool(x1,y1)
    % a*x^2+b*x-d-exp(x*c)
    % cftool(double(log(x_scaled(mask))),double(f_x_scaled(mask)))
    % 
    % x1=double(log10(x_scaled(mask)));
    % y1=double(f_x_scaled(mask));
    % T = table(x1(:), y1(:), 'VariableNames', {'x', 'y'});
    % writetable(T, ['data_' num2str(k) '.csv']);

    sh_str = sprintf('%.2e', double(E0(k)));  % display only, keep double for label
    legendEntries{end+1} = ['$E_0 = ' sh_str '$'];
    k = k + 1;
end

ylim([1e-15 1e10 ])
xlim([1e-10, 1e3])

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';

xlabel(ax, '$x$', ...
       'Interpreter', 'latex');
ylabel(ax, '$ f(x) $', ...
       'Interpreter', 'latex');
xticks([1e-10 1e-7 1e-4 1e-1 1e2])
% ylabel(ax, '$ \ln{f(E)} + \ln {\langle E \rangle} - \frac{E}{\langle E \rangle} $', ...
%        'Interpreter', 'latex');
legend(legendEntries, 'Interpreter', 'latex', ...
    'Location', 'northeast','FontSize',14);

ax3 = axes('Position',[0.21 0.25 0.28 0.28]);
h=plot(time,meanEnerg);
xlim([1 1e4])

h.LineWidth=1;
ax=gca;
% ax.XScale='log';
ax.YScale='log';
ax.FontSize=12;

% Set custom ticks at specific values
xticks([ 1000, 10000])
yticks([10^0,10^10, 10^20])
% Set custom tick labels (LaTeX-formatted)
xticklabels({'$ 10^3$', '$10^4$'})

% Enable LaTeX interpreter
ax = gca;
ax.TickLabelInterpreter = 'latex';

xlabel('$t$',Interpreter='latex')
ylabel('$\langle E_k \rangle$',Interpreter='latex')
t=text(3.5e3,1e10,'$\langle E_k \rangle \propto \exp(t)$',Interpreter='latex');
t.FontSize=12;

exportgraphics(fig1, 'Ultra_relativistic_self_similarities_V2.png', 'Resolution', 300)