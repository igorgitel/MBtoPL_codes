function [v_b0,f_sim0,f_hist] = histlog(data,edges)
    %[f_sim0,edges] = histcounts(data,edges,'Normalization','countdensity');
    [f_sim0,edges] = histcounts(data,edges);
    f_hist=f_sim0;
    Le=log10(edges);
    centers=zeros(1,length(edges)-1);
    for i=1:length(edges)-1
        centers(i)=(Le(i)+Le(i+1))/2;
        if f_sim0(i)<50
            f_sim0(i)=0;
        else           
            f_sim0(i)=f_sim0(i)/(edges(i+1)-edges(i));
        end
    end
    v_b0=10.^(centers);
end