function [e1, p1,e2, p2] = RelCol(e1, p1,e2, p2,n0)
%collision is function that colide two particles
% Center mass momentum and energy

V=(p1+p2)/(e1+e2);

if V==p1/e1
    error("V==p1/e1")
else if V==p2/e2
    error("V==p2/e2")
end

V=double(V);

if norm(V)>1
    error('V larger then 1')
end

if prod(~isreal(V))==1
    error('V isnt real')
end

[e10, p10]=LorentzTransform(e1,p1,V);

if ~isreal(e10)==1
    error('e10 isnt real')
end

[e20, p20]=LorentzTransform(e2,p2,V);

if ~isreal(e20)==1
    error('e10 isnt real')
end

Np10=norm(p10);
Np20=norm(p20);

if ~isreal(Np10)==1
    error('Np10 isnt real')
end


if ~isreal(Np20)==1
    error('Np10 isnt real')
end
%  p10+p20;
% if (p10+p20)>=1e-12
%     error('p10+p20~=0')
% end
% assign the directions to the particles velocity
% in center of mass frame after the collision (0 - c.m. quantity)
% n0=randdir;
p01after=Np10*n0;
p02after=-Np20*n0;
e102=sqrt(p01after*p01after'+1);
e202=sqrt(p02after*p02after'+1);

% Np10-Np20
% Transformation back to lab frame
[e1,pL1a]=LorentzTransform(e10,p01after,-V);
if ~isreal(e1)==1
    error('e1 isnt real')
end
[e2,pL2a]=LorentzTransform(e20,p02after,-V);
if ~isreal(e2)==1
    error('e2 isnt real')
end
% if abs(e1-sqrt(dot(pL1a,pL1a)+1))>1e-10
%     error('e1 no equal to E(p1)')
% end
% if abs(e2-sqrt(dot(pL2a,pL2a)+1))>1e-10
%     error('e2 no equal to E(p2)')
% end
%% assign
%     pLab1after1=pLab1after(1); pLab1after2=pLab1after(2); pLab1after3=pLab1after(3);
%     pLab2after1=pLab2after(1); pLab2after2=pLab2after(2); pLab2after3=pLab2after(3);
%

p1=pL1a;
p2=pL2a;
end