function [bins,fsim10,fsim50,fsim100,fsim,f_hist] = low_err_histlog(data,edges)
    %[f_sim0,edges] = histcounts(data,edges,'Normalization','countdensity');
    [f_sim0,edges] = histcounts(data,edges);
    f_hist=f_sim0;
    Le=log10(edges);
    centers=zeros(1,length(edges)-1);

    for i=1:length(edges)-1
        centers(i)=(Le(i)+Le(i+1))/2;

        fsim(i)=f_sim0(i)/(edges(i+1)-edges(i));

        if f_sim0(i)<10
            fsim10(i)=0;
            fsim50(i)=0;
            fsim100(i)=0;
        else 
            fsim10(i)=f_sim0(i)/(edges(i+1)-edges(i));
            if f_sim0(i)<50
                fsim50(i)=0;
                fsim100(i)=0;
            else
                fsim50(i)=f_sim0(i)/(edges(i+1)-edges(i));
                if f_sim0(i)<100
                    fsim100(i)=0;
                else
                    fsim100(i)=f_sim0(i)/(edges(i+1)-edges(i));
                end
            end
        end
    end
    bins=10.^(centers);
%     trapz(bins,fsim)
%     trapz(bins,fsim10)
%     trapz(bins,fsim50)
%     trapz(bins,fsim100)
end