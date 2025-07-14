function [pLab1a,Ea1] = Inf_Rel_Col(p1,e1,beta2)
%collision is function that colide two particles 
%     
%% Momenrum definition
% p1=[px1 py1 pz1];
%% Rest frame of the heavy particle momentum and energy
    [e10, p10]=LorentzTransform(e1,p1,-beta2);
%% Assign the directions to the particles velocity
    % in center of mass frame after the collision (0 - c.m. quantity)
    e10after=norm(p10);
    p10after=e10after*randdir;
%% Transformation back to lab frame
    [Ea1,pLab1a]=LorentzTransform(e10after,p10after,beta2);
end


