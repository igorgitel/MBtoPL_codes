function n = randdir_mat(N)
% randomly choose the diration
    phir=2*pi*rand(N,1);
    cosThetar=2*rand(N,1)-1;
    sinThetar=sqrt(1-cosThetar.^2);
% randomly chosen direction
    n(:,1)=sinThetar.*cos(phir); 
    n(:,2)=sinThetar.*sin(phir);
    n(:,3)=cosThetar;
end