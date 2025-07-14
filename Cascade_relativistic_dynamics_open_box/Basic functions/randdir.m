function n = randdir()
%randdir randomly choose direction. 
phi=2*pi*rand;
cosTheta=2.*rand-1;
sinTheta=sqrt(1-cosTheta.^2);
%% assign 
n(1)=sinTheta.*cos(phi);
n(2)=sinTheta.*sin(phi);
n(3)=cosTheta;
end