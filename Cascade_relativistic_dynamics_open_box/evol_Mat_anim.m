function evol_Mat_anim(data,str,folder,ic,icb)
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
nImages = length(data.time)-1;
for idx = 1:nImages

    plot(data.bins,data.fsim10(idx,:),'o-','LineWidth',2);
%     hold on
% %     plot(data.bins,data.Teo,'-','LineWidth',2);
%     hold off
    % title
    title(['Time step ' num2str(data.time(idx),'%.0f')], 'Interpreter', 'latex', 'FontSize', 16)

    % labels
    xlabel( 'The kinetic energy $E$', 'Interpreter', 'latex', 'FontSize',18);
    ylabel( '$f(E)$', 'Interpreter', 'latex', 'FontSize',18);

    % Limits
    ylim([1e-5 1e30])
    xl=[1e-10 1e8];
    xlim(xl)

    %tickets
    xticks(10.^(-16:3:16))
    yticks(10.^(-16:3:16))

    ax = gca;
    ax.FontSize = 12;
    ax.XScale='log';
    ax.YScale='log';
    %% 
    xline(ic.ke,'-.',{'Initial', 'kinetic', 'energy'})
    xline(icb.ke,'-.',{'Final', 'kinetic', 'energy'})

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
disp('Animatetion Done')
end