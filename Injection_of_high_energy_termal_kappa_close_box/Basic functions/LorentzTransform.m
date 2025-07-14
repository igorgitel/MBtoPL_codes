function [Ea,Pa] = LorentzTransform(E,p,V)
%Lorentz transformation of (E,p) to frame thats move with Velocity V.
%
% digits(1000);
% E=vpa(E);
% p=vpa(p);
% V=vpa(V);

b=norm(V);
if ~isreal(b) || b>=1
    error('beta isnt real')
end

g=1/sqrt(1-b^2);

if ~isreal(g)
    error('gamma isnt real')
end

% 
p4=[E p]';
gm=(g-1)/(V*V');
L = [g -g*V(1) -g*V(2) -g*V(3);...
    -g*V(1) 1+gm*V(1)^2 gm*V(1)*V(2) gm*V(1)*V(3);...
    -g*V(2) gm*V(2)*V(1) 1+gm*V(2)^2 gm*V(2)*V(3);...
    -g*V(3) gm*V(3)*V(1) gm*V(3)*V(2)  1+gm*V(3)^2];

pT=L*p4;

%% The transformation of E and the vector p
% if g<10^10
% n=V/b;
% Ea=g*(E-b*dot(n,p));
% Pa=p+(g-1)*dot(n,p)*n-g*E*b*n;
% 
% abs(Ea^2-(dot(Pa,Pa)+1))
% abs(pT(1)-Ea)

% if abs(Ea^2-(dot(Pa,Pa)+1))>1e-10
%     error('e(p) not equal e')
% end

% 
% Efp=sqrt(dot(Pa,Pa)+1);
% dif=Efp-Ea;

Ea=double(pT(1));
Pa=double(pT(2:4))';

Efp=sqrt(dot(Pa,Pa)+1);
% dif=Efp-Ea

% if abs(Ea^2-Efp)>1e-10
%     error('E dosnt equal E(pT)')
% end
end