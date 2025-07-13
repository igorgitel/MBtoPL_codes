function [vx, vy, vz] = ChooseRandomV(v0,N)
%Choose random diraction of the speed
%   v0 - the module of the speed 
%   N - number of particle
Ranprep=rand(N,2);
phi(:)=2*pi*Ranprep(:,1);
cosTheta(:)=2.*Ranprep(:,2)-1;
sinTheta(:)=sqrt(1-cosTheta(:).^2);
%% assign 
vx(:)=v0.*sinTheta(:).*cos(phi(:));
vy(:)=v0.*sinTheta(:).*sin(phi(:));
vz(:)=v0.*cosTheta(:);
end

