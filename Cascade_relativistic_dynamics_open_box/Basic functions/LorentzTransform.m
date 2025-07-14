function [E,pT1] = LorentzTransform(E,p,V)
% Lorentz transformation of (E,p) to frame thats move with Velocity V. And
% get the four-vector (E,pT1). 
% 
beta=norm(V); %beta=v/c
gam=1/sqrt(1-beta^2);%the gamma factor
% gam=1/sqrt(1-beta^2);
if beta==0
    n=[0 0 0]; % the direction in 3D
else
    n=V/beta; % the direction in 3D
end
p4=[E p]; % the 4 vector before the transformation 

L=[gam -gam*beta*n(1) -gam*beta*n(2) -gam*beta*n(3);...
    -gam*beta*n(1) 1+(gam-1)*((n(1))^2) (gam-1)*n(1)*n(2) (gam-1)*n(1)*n(3);...
    -gam*beta*n(2) (gam-1)*n(2)*n(1) 1+(gam-1)*((n(2))^2) (gam-1)*n(2)*n(3);...
    -gam*beta*n(3) (gam-1)*n(3)*n(1) (gam-1)*n(3)*n(2)  1+(gam-1)*((n(3))^2)];

% the lorentz matrix

pT=L*(p4.'); % multiplication of matrix L by the column 4-vector. (' stand for transpose)

E=double(pT(1)); % assign the energy 
pT1=double(pT(2:4))';% assign the momentUM

if isnan(E)
    error('error E is nun in Lorentz')
end

end

