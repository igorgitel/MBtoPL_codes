function f = frequancy(v1,v2,cs_p)
% frequancy calculation
v=norm(v1-v2);
f=sigma(v,cs_p)*v;
end