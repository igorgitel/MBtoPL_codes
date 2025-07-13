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
% Np10-Np20
% Transformation back to lab frame
[e1,p1]=LorentzTransform(e10,p01after,-V);
if ~isreal(e1)==1
    error('e1 isnt real')
end
[e2,p2]=LorentzTransform(e20,p02after,-V);
if ~isreal(e2)==1
    error('e2 isnt real')
end

end