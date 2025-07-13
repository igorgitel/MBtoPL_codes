function [e1a,e2a]...
            = collision(v1,v2,m1,m2,n0)
%collision is functionthat colide two particles 

%% v Before
   v=norm(v1-v2);
%% assign the directions to the particles velocity
    % in center of mass frame
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

