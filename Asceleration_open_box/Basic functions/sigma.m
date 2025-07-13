function sigma = sigma(v,cs_p)
% cross section calculation
if v<0.1
    sigma=(0.1).^(cs_p);
else
    sigma=v^(cs_p);
end