function n = randdir()
% randomly choose the diration
    phir=2*pi*rand;
    cosThetar=2*rand-1;
    sinThetar=sqrt(1-cosThetar^2);
% randomly chosen direction
    n(1)=sinThetar.*cos(phir); 
    n(2)=sinThetar.*sin(phir);
    n(3)=cosThetar;
end