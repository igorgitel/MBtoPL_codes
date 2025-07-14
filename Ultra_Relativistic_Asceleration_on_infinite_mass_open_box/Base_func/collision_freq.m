function f_now = collision_freq(p1,e1,beta2)

% Calculation the colision frequamcy
v1=p1./e1;
v2=beta2;

    % Theoretical calculation of collision frequancy.
    Cr=cross(v1,v2);
    Res=v1-v2;
%     v_rel=sqrt(dot(Res,Res)-dot(Cr,Cr))/(1-dot(v1,v2));
    sigma=1;
    f_now=sigma*sqrt(dot(Res,Res)-dot(Cr,Cr));
end
