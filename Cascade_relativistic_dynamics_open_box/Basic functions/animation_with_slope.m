function animation_with_slope(s,str,ic,ts_jump,y0_line_value)

FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

f=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);
% f.Position([1 2])=[100 100];
% f.Position([3 4])=[1080 720];

axis tight manual
set(gca,'nextplot','replacechildren');
%%
% load('m_1.0_MB_Velocities_N_1e+05.mat','vx','vy','vz');
% energ=(vx.^2+vy.^2+vz.^2)/2;
% edges=10.^(-16:0.1:16);
% [e_b0,f_sim0] = histlog(energ,edges);
% f_sim0=f_sim0./length(energ);
% clear vx vy vz energ edges
%%
% Create a video writer object for the output video file and open the object
% for writing.
%filename=['Number_of_par' N1_str 'coll_num_' Total_coll_num_str '_p0_' p0_str '.mp4'];

filename=[fullfile('Videos_02',str) '.mp4'];
v = VideoWriter(filename,'MPEG-4');
open(v);
%%
% Generate a set of frames, get the frame from the figure, and then write each
% frame to the file.
nImages = length(s);
k=0;
for idx = 1:ts_jump:nImages

    k=k+1;
    sb1=subplot(4,1,[1 3]);
    sb1.Position=[0.1300    0.34    0.7750    0.5959];
    plot(s(idx).ebin,s(idx).fsim,'-o','LineWidth',2);

    % labels
    ylabel( '$F(E_k)$', 'Interpreter', 'latex');

    % Limits
    yl=[1e-7 1e1];
    ylim(yl)
    xl=[1e2 1e5];
    xlim(xl)

    %tickets
    xticks(10.^((log10(xl(1))+1):1:(log10(xl(2))-1)))
    yticks(10.^((log10(yl(1))+1):1:(log10(yl(2)))))

    % axis properties
    grid off
    ax = gca;
    ax.FontSize = 14;
    ax.XScale='log';
    ax.YScale='log';

    % slope
    sb2=subplot(4,1,4);
    sb2.Position=[0.1300    0.1100    0.7750    0.23];

    h=3;
    i=1;
    while i<=h
        y(i)=log10(s(idx).fsim(i+h))-log10(s(idx).fsim(i));
        y(i)=y(i)/(log10(s(idx).ebin(i+h))-log10(s(idx).ebin(i)));
        x(i)=s(idx).ebin(i);
        i=i+1;
    end
    for i=(h+1):(length(s(idx).ebin)-(h+1))
        y(i)=log10(s(idx).fsim(i+h))-log10(s(idx).fsim(i-h));
        y(i)=y(i)/(log10(s(idx).ebin(i+h))-log10(s(idx).ebin(i-h)));
        x(i)=s(idx).ebin(i);
    end
    while i<=length(s(idx).ebin)
        y(i)=log10(s(idx).fsim(i))-log10(s(idx).fsim(i-h));
        y(i)=y(i)/(log10(s(idx).ebin(i))-log10(s(idx).ebin(i-h)));
        x(i)=s(idx).ebin(i);
        i=i+1;
    end

    p=plot(x,y,'-','LineWidth',2,'MarkerSize',2);
    p.LineWidth=2;
    xlim(xl)
    ylim([-2.7 1])
    xticks(10.^((log10(xl(1))):2:(log10(xl(2)))));
    yticks([-3 -2.5 -2 -1.5 -1 -0.5 0])
    xlabel('$E_k$',Interpreter='latex');

    ylabel('PL index')
    % yticks(-6:0.5:-2);
    
    y_line=yline(y0_line_value);
    y_line.LineWidth=4;

    % grid properties
    grid off
    ax = gca;
    ax.FontSize = 12;
    ax.XScale='log';


    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
display('Animatetion Done')
end