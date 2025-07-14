function [s] =...
    TwotypeColision(N,ic,edges,time_s)
%% Choose the direction
[px,py,pz,energ] = ChooseRandomRelV(ic.p01,N);
%% Calculation
coldone=0; % number of collision done so far
step=0;
time=0;
f_max=0;
idx=0;
maxplotnum=5e3;
Energy_mean=ic.E1;
temp1=0;
s=struct('EnergyBE1',cell([1 maxplotnum]),...
    'DFE1',cell([1 maxplotnum]),...
    'coldone',cell([1 maxplotnum]),...
    'coltry',cell([1 maxplotnum]),...
    'Energ_mean',cell([1 maxplotnum]));
% while time<time_s.maxtime
while Energy_mean<1e35
    %% Collision
    pID1=round(N*rand+0.5);

    p1=[px(pID1) py(pID1) pz(pID1)]; e1=energ(pID1);
    beta2=ic.beta2*randdir;
%     beta2=ic.beta2*[0 1 0];
    f_now =collision_freq(p1,energ(pID1),beta2);
    if f_max==0
        f_max=f_now;
    end
    time=time+(1/N)*(1/f_max);
    f_max=max([f_max f_now]);
    if f_now>f_max*rand
        % Kinematics of the collision
        [p1,e1] = Inf_Rel_Col(p1,e1,beta2);
        
        px(pID1)=p1(1); py(pID1)=p1(2); pz(pID1)=p1(3);
        coldone=coldone+1;
        energ(pID1)=e1; 
    end
     if time>step
        step=step+time_s.stepinc; 
        idx=idx+1;
        [s(idx).EnergyBE1,s(idx).DFE1] = ...
            histlog(energ,edges.type1);
      
        s(idx).DFE1=s(idx).DFE1./N;
        col_this_time_step=coldone-temp1;
        temp1=coldone;
        s(idx).coldone=coldone;
        s(idx).time=time;
        Energy_mean=mean(energ(1:N));
        s(idx).Energ_mean=Energy_mean;
        if mod(time,1e1)<1e-1
            disp({'time step' 'ts col' 'Mean energy';...
                time col_this_time_step  Energy_mean})
        end
     end
end
sf=s(1:floor(time));
clear s 
s=sf;
end

