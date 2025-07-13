function f_teo = Teo_kin_energy(edges,Mass,T_MJ)
% Calculation the therotical equlibrium distribution function in kinetic
% energy
    % The distribution in energy
    fMJE = @(x,m) 1/(m.^2*T_MJ*besselk(2,m/T_MJ)).*x.*sqrt(x.^2-m.^2)...
        .*exp(-(x./T_MJ));
    % The distribution in kinetic energy
    fMJ_kE=@(x,m,T) 1/(m.^2.*T.*besselk(2,m/T,1)).*(x+m)...
        .*(sqrt((x+m).^2-m.^2)).*exp(-x/T);
    % The distribution in momentum devided by the mass
    fMJmu = @(x,m) m/(T_MJ.*besselk(2,m/T_MJ))...
        .*x.^2.*exp(-(m./T_MJ)*sqrt(1+x.^2));

   f_teo=ones(1,length(edges)-1);
%    centers=zeros(1,length(edges)-1);
%     for i=1:length(edges)-1
%         centers(i)=(Le(i)+Le(i+1))/2;
%     end
    for i=1:(length(edges)-1)
        q1=integral(@(x) fMJ_kE(x,Mass,T_MJ),edges(i),edges(i+1));
        if isreal(q1)
            f_teo(i) = q1./(edges(i+1)-edges(i));
        else
            f_teo(i)=NaN;
        end
    end
%     plot(centers,f_teo)
end