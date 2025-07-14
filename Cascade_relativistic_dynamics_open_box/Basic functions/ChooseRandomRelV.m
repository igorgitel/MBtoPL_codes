function [px,py,pz,e] = ChooseRandomRelV(p0,N,m)
%Choose random diraction of the speed
%   v0 - the module of the speed 
%   N - number of particle
Ranprep=rand(N,2);
phi(:)=2*pi*Ranprep(:,1);
cosTheta(:)=2.*Ranprep(:,2)-1;
sinTheta(:)=sqrt(1-cosTheta(:).^2);
%% assign 
px(:)=p0.*sinTheta(:).*cos(phi(:));
py(:)=p0.*sinTheta(:).*sin(phi(:));
pz(:)=p0.*cosTheta(:);
%% reverse diractions half
% px0=-px; px=[px px0];
% py0=-py; py=[py py0];
% pz0=-pz; pz=[pz pz0];

pnorma(:)=sqrt(px(:).^2+py(:).^2+pz(:).^2);
e(:)=sqrt(pnorma(:).^2+m^2);
end

