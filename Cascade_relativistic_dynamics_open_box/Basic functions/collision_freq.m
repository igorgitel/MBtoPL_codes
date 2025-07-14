function f_now = collision_freq(p1,e1,p2,e2,f_param,cs_power)
% Calculation the colision frequancy


    p41co=[e1 p1];
    p42contra=[e2 -p2]';
    % The velocityies of the particles 
%     v1=p1./e1;
    v2=p2./e2;
    
    % the relative velocity for cross section
    [E_rel,~] = LorentzTransform(e1,p1,v2);
%     prel=norm(p_rel);
%     vrel=norm(p_rel)/E_rel;
    % Theoretical calculation of collision frequancy.
%     Cr=cross(v1,v2);
%     Res=v1-v2;
if E_rel<0.1
    sigma=(0.1)^(cs_power);
else
    sigma=(E_rel)^(cs_power);
end
    % calculation
%      f_now=f_param*sigma*sqrt(dot(Res,Res)-dot(Cr,Cr));   
    f_now=f_param*sigma*(sqrt((p41co*p42contra)^2-1))/(e1*e2);
    if ~isreal(f_now)
        error('frequancy isnt real')
    end
end