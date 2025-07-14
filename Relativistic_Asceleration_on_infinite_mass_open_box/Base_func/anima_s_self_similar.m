function anima_s_self_similar(s,p0,str,N,edges,T,m1,m2)
f=figure

f.Position([1 2])=[100 100];
f.Position([3 4])=[1280 720];

axis tight manual 
set(gca,'nextplot','replacechildren'); 
%% 
% Create a video writer object for the output video file and open the object 
% for writing.

Total_coll_num_max=max(j)*N(1);
Total_coll_num_str=num2str(Total_coll_num_max,'%1.0e');
p0_str=num2str(p0.m1,'%1.3f');
N1_str=num2str(N(1),'%1.0e');
%filename=['Number_of_par' N1_str 'coll_num_' Total_coll_num_str '_p0_' p0_str '.mp4'];
filename=['self_similar_' str '_p0_' p0_str  '.mp4'];
v = VideoWriter(filename,'MPEG-4');
open(v);
%% 
% Generate a set of frames, get the frame from the figure, and then write each 
% frame to the file.
nImages = length(s)-1;
 k=0;
for idx = 1:nImages
       k=k+1;
       MeanE=trapz(s(idx).EnergyBE1,s(idx).EnergyBE1.*s(idx).DFE1);
%         MeanE=s(idx).Energ_mean;
        MaxF=max(s(idx).DFE1);
       plot(s(idx).EnergyBE1./MeanE,s(idx).DFE1./MaxF,'o-','LineWidth',2);
%      plot(s(end).EnergyBE1,s(end).Teo1,'-','LineWidth',2,'Color','k');
       hold off
      
        % title
        title({['$\gamma_1 \beta_1 =$ ' num2str(p0.gamma_beta,'%1.1e')]...
            ,['Number of collisions performed '...
            num2str(s(idx).coldone/N,'%1.1e')],...
... ['Time is ' num2str(s(idx).time,'%0.1f') ],...
              [' The mass ratio is ' num2str(m1/m2,'%1.0e')]},...
              'Interpreter', 'latex', 'FontSize', 16) 
        
        % labels
        xlabel( '$\frac{\mathcal{E}}{\mathcal{E}_{mean}}$', 'Interpreter', 'latex', 'FontSize',22);
        ylabel( '$\frac{f}{f_{max}}$', 'Interpreter', 'latex', 'FontSize',22);

        
        %legends
        leg={['Simulated DF for m1=' num2str(m1)]};
        lg=legend(leg,'Location','southwest','NumColumns',1);
        title(lg,['N = ' num2str(N)]);
        
        
        % Limits
        ylim([1e-5 10])
        xlim([1e-4 1e4])
        
        % grid properties
            grid on
            grid minor
            ax = gca;
            ax.FontSize = 16;
           ax.XScale='log';
           ax.YScale='log';
            
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);
display('Animatetion Done')
end