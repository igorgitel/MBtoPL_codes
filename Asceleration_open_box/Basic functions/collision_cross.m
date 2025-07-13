function [e1a,e2a]...
            = collision_cross(v1,v2,m1,m2,n0)
%collision is functionthat colide two particles 

%% v Before
   v=norm(v1-v2);
%% assign the directions to the particles velocity
    % in center of mass frame
%     max_dsigma=0.41;
% dsigma=@(x) 3/8*(1+cos(x)^2)*sin(x);
% 
% num_trys=0;
% chek=false;
% while chek==false
%     n0=randdir;
%     theta=n0*p1'/norm(p1);
%     if dsigma(theta)<max_dsigma*rand
%         chek=true;
%     end
%     num_trys=num_trys+1;
% end
%     n0=randdir;
    vcm1=(m2/(m1+m2))*v.*n0;
    vcm2=-(m1/(m1+m2))*v.*n0;
%% Center of mass velocity
    V=(m1*v1+m2*v2)/(m1+m2);
%% Transformation back to lab frame
    v1=vcm1+V;
    v2=vcm2+V; 
    %% The energy of the particles. 
    e1a=m1*(v1*v1')/2;
    e2a=m2*(v2*v2')/2;
end

