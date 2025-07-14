function anima_s_2(s,str,folder)
figure;

axis tight manual
set(gca,'nextplot','replacechildren');
%%
% Create a video writer object for the output video file and open the object
% for writing.
filename=[fullfile(folder,str) '.mp4'];
v = VideoWriter(filename,'MPEG-4');
open(v);
%%
% Generate a set of frames, get the frame from the figure, and then write each
% frame to the file.
nImages = length(s)-1;
for idx = 1:nImages
    subplot(4,1,[1 3])
    plot(s(idx).bins,s(idx).f_sim,'o-','LineWidth',2);
    hold on
    plot(s(idx).bins,s(idx).Teo,'-','LineWidth',2);
    hold off
    % title
    title(['Time step ' num2str(s(idx).time)], 'Interpreter', 'latex', 'FontSize', 16)

    % labels
    xlabel( 'The kinetic energy $E$', 'Interpreter', 'latex', 'FontSize',18);
    ylabel( '$f(E)={1 \over N}{dN \over dE}$', 'Interpreter', 'latex', 'FontSize',18);

    % Limits
    ylim([1e-5 1e10])
    xl=[1e-15 1e-2];
    xlim(xl)

    %tickets
    xticks(10.^(-16:2:16))
    yticks(10.^(-16:2:16))

    % grid properties
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12;
    ax.XScale='log';
    ax.YScale='log';

    x=s(idx).bins; y=s(idx).f_sim;
    subplot(4,1,4)
    h=3;
    for i=(h+1):(length(y))
        dy(i)=log10(y(i))-log10(y(i-h));
        dy(i)=dy(i)/(log10(x(i))-log10(x(i-h)));
    end
    plot(s(idx).bins,dy,'o-','LineWidth',2);

    xlabel( 'The kinetic energy $E$', 'Interpreter', 'latex', 'FontSize',18);
    ylabel('Derivative')
    % Limits
    ylim([-1.5 0])
    xlim(xl)

    % grid properties
    grid on
    grid minor
    ax = gca;
    ax.FontSize = 12;
    ax.XScale='log';
    % ax.YScale='log';
    
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
disp('Animatetion Done')
end