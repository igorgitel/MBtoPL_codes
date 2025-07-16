function anima_s(data,str,folder)
FigureSize = [10 10 21 13];
DefaultFontSizeForFigure=14;

fig1=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);
axis tight manual
set(gca,'nextplot','replacechildren');
%%
% Create a video writer object for the output video file and open the object
% for writing.
% filename=[fullfile(folder.one,folder.two,str) '.mp4'];
filename=[fullfile(folder,str) '.mp4'];
v = VideoWriter(filename,'MPEG-4');
open(v);
%%
% Generate a set of frames, get the frame from the figure, and then write each
% frame to the file.
nImages = length(data)-1;
for idx = 1:nImages


    sb1=subplot(4,1,[1 3]);
    sb1.Position=[0.1300    0.34    0.7750    0.5959];
    p=plot(data(idx).bins,data(idx).f_sim,"-o",LineWidth=2);
    p.MarkerSize=2;
    hold on
    plot(data(idx).bins,data(idx).Teo,"-",LineWidth=1);
    hold off
    ax=gca;
    ax.XScale='log';
    ax.YScale='log';
%     xl=[1e-6 1e2];
%     yl=[1e-4 1e6];
    xl=[1e-16 1e2];
    yl=[1e-16 1e16];
    xlim(xl)
    ylim(yl)

    xticks(10.^[(log10(xl(1))):(log10(xl(2)))])
    xticklabels([])
    yticks(10.^[(log10(yl(1))+1):1:(log10(yl(2))-1)])

    xlabel('$Mass$',Interpreter='latex')
    ylabel('$f(m)$',Interpreter='latex')
    grid off
%     title(['TS ' num2str(data(idx).time,'%.0f') ]);

    x=data(idx).bins; y=data(idx).f_sim;

    %% time step
    annotation('textbox', [0.6, 0.8, 0.2, 0.1], 'String',...
        ['Time t = ' num2str(data(idx).time,'%.2f')], 'FitBoxToText', 'on',...
        'FontSize', 16, 'FontWeight', 'bold', ...
        'BackgroundColor', [0.8, 0.8, 0.8], 'Interpreter', 'latex', ...
        'EdgeColor', 'none');

%     t = text(xl(2)/100,yl(2)/10,['Time t = ' num2str(data(idx).time,'%.2f')]);
%     t.FontSize=14;
    %%
    sb2=subplot(4,1,4);
    sb2.Position=[0.1300    0.1100    0.7750    0.23];
    % sb2.OuterPosition=[0    0.02   1.0000    0.335];
    h=3;
    for i=(h+1):(length(y))
        dy(i)=log10(y(i))-log10(y(i-h));
        dy(i)=dy(i)/(log10(x(i))-log10(x(i-h)));
    end
    h=plot(x,dy,'o-','LineWidth',1.5);
    h.MarkerSize=3;
    ax=gca;
    ax.XScale='log';
    ylim([-4 1])
    xlim(xl)
    % y0=mean(dy(x>1e1 & x<5e2));
    % y_std=std(dy(x>1e1 & x<5e2));
    % yline(y0,'-',[' ' num2str(y0,'%.2f') '$\pm$' num2str(y_std,'%.2f')],LineWidth=1,Interpreter='latex');
    xlabel('$Mass$',Interpreter='latex')
    ylabel('PL index',Interpreter='latex')
    xticks(10.^[(log10(xl(1))):(log10(xl(2)))])


    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
disp('Animatetion Done')
end