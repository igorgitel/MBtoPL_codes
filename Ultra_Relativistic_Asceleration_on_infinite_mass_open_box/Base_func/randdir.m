function n0 = randdir
% return a random 3 dimantional diraction vector, with unit magnitude.
    phir=2*pi*rand;
    cosThetar=2*rand-1;
    sinThetar=sqrt(1-cosThetar^2);
% randomly chosen direction
    n0=[sinThetar.*cos(phir) sinThetar.*sin(phir) cosThetar];
end
