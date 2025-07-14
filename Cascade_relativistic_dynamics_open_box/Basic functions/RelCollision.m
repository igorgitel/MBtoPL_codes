function [pLab1after,Eafter1,...
            pLab2after,Eafter2]...
            = RelCollision(p1,e1,...
                           p2,e2,m1,m2)
%collision is function that colide two particles 
%     
%% Momenrum definition
% p1=[px1 py1 pz1];
% p2=[px2 py2 pz2];
%% Center mass momentum and energy
    V=(p1+p2)/(e1+e2);

    [e10, p10]=LorentzTransform(e1,p1,V);
    [e20, p20]=LorentzTransform(e2,p2,V);
% 
    Np10=sqrt(p10*p10');
    Np20=sqrt(p20*p20');

%% randomly choose the diration of motion in CM
    phir=2*pi*rand;
    cosThetar=2*rand-1;
    sinThetar=sqrt(1-cosThetar^2);
% randomly chosen direction
    n0=[sinThetar.*cos(phir) sinThetar.*sin(phir) cosThetar];
%     norm(n0);

%% assign the directions to the particles velocity
    % in center of mass frame after the collision (0 - c.m. quantity)
    p01after=Np10*n0;
    p02after=-Np10*n0;
    e10after=sqrt(norm(p01after)^2+m1^2);
    e20after=sqrt(norm(p02after)^2+m2^2);
%% Transformation back to lab frame
    [Eafter1,pLab1after]=LorentzTransform(e10after,p01after,-V);
    [Eafter2,pLab2after]=LorentzTransform(e20after,p02after,-V); 
end

