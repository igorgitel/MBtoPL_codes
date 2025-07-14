function anima_s(s,str,folder)
f=figure;

f.Position([1 2])=[100 100];
f.Position([3 4])=[1080 720];

axis tight manual 
set(gca,'nextplot','replacechildren'); 
%% 
% Create a video writer object for the output video file and open the object 
% for writing.
%filename=['Number_of_par' N1_str 'coll_num_' Total_coll_num_str '_p0_' p0_str '.mp4'];
filename=['zzh_' str '.mp4'];
v = VideoWriter(fullfile(folder,filename),'MPEG-4');
open(v);
%% 
% Generate a set of frames, get the frame from the figure, and then write each 
% frame to the file.
nImages = length(s)-1;
 k=0;
 for idx = 1:nImages
       k=k+1;
        plot(s(idx).ebin,s(idx).fsim,'o-','LineWidth',2);
        % title
        title(['Time step ' num2str(idx) 'particle number is ' ...
            num2str(s(idx).parNum)], 'Interpreter', 'latex', 'FontSize', 16) 
        
        % labels
        xlabel( 'The kinetic energy $E$', 'Interpreter', 'latex', 'FontSize',18);
%         ylabel( '$f(E)={1 \over N}{dN \over dE}$', 'Interpreter', 'latex', 'FontSize',18);
        ylabel( '$f(E)$', 'Interpreter', 'latex', 'FontSize',18);
        
        % Limits
        ylim([1e-12 1e20])
        xlim([1e-15 1e0])

        %tickets 
        xticks(10.^(-16:3:16))
        yticks(10.^(-16:3:16))

        % grid properties
            grid on
            grid minor
            ax = gca;
            ax.FontSize = 32;
           ax.XScale='log';
           ax.YScale='log';
            
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);
display('Animatetion Done')
end