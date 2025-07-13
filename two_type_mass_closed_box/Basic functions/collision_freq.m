function f_now = collision_freq(p1,e1,p2,e2,alpha)

% Calculation the colision frequamcy
v1=p1./e1;
v2=p2./e2;

Cr=cross(v1,v2);
Res=v1-v2;
[e_rel,p_rel] = LorentzTransform(e1,p1,v2);
v_rel=abs(p_rel)/e_rel;
% sig=sigma(v_rel,alpha);
sig=1;
f_now=sig*sqrt(dot(Res,Res)-dot(Cr,Cr));
if ~isreal(prod(f_now))
    error('f isnt real');
end
end