clear

folder='saves';
gb=1e-3;
fparam=1e1;
m_ratio=1836;
max_time=2e4
N=1e6;
str=['Mass_ratio_' num2str(m_ratio)...
    '_gb_' num2str(gb,'%.0e') ...
    '_f_param_' num2str(fparam,'%.0e')...
    '_Nop_' num2str(N) ...
    '_max_time_' num2str(max_time) '.mat'];
% str=['Mass_ratio_' num2str(m_ratio)...
%     '_gb_' num2str(gb,'%.0e') ...
%     '_f_param_' num2str(fparam,'%.0e')...
%     '_Nop_' num2str(N) '.mat'];
load(fullfile(folder,str))

length(data.f_sim_1) 
ind=[20 550 1000 3000 5000 1e4 2e4];
%%
FigureSize = [0 0 21 13];
DefaultFontSizeForFigure=14;

fig=figure('Units','centimeters','Position',FigureSize,...
    'DefaultAxesFontSize',DefaultFontSizeForFigure);

for i_relevant=1:length(ind)

    for i=1:length(ind)
        if i==i_relevant;
            x=data.bins;
            y1=data.f_sim_1(ind(i),:);
            y1=y1/trapz(x,y1);
            plot(x,y1,'r:.',LineWidth=3)
            hold on

            y2=data.f_sim_2(ind(i),:);
            y2=y2/trapz(x,y2);
            plot(x,y2,'b-.',LineWidth=3);

            if i==length(ind)
                title({['$t_{eq}=' num2str(data.time(ind(i)),'%.0f') ',\,\,N_{cpp}=' num2str(data.coll_num_Mm(ind(i))/(2*N),'%.1f') '$']},...
                    Interpreter='latex',FontSize=26,FontWeight='bold')
            else
                 title({['$t=' num2str(data.time(ind(i)),'%.0f') ',\,\,N_{cpp}=' num2str(data.coll_num_Mm(ind(i))/(2*N),'%.1f') '$']},...
                    Interpreter='latex',FontSize=26,FontWeight='bold')
            end
        end
    end
    plot(data.bins,data.Teo1,'k-')
    % plot(data.bins,data.Teo2,'y-')

    ax=gca;
    ax.XScale='log';
    ax.YScale='log';
    ax.FontSize=16;
    xlabel('$E_k$',Interpreter='latex')
    ylabel('$f(E_k)$',Interpreter='latex')

    ylim([1e-5 1e10])
    xlim([1e-8 1e-1])

    grid off

    h=annotation("textarrow", [0.2293 0.1902], [0.2617 0.4053], "String","MB");
    h.FontSize=16;
    h.Interpreter='latex';

    hold off
    exportgraphics(fig,['Two_types_of_particle_e_' num2str(i_relevant) '.png'],'Resolution',300)
end

%% 
%
