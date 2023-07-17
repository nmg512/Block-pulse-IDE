function F = threeratesol(L, j, m, gj, gj1, gj2, r)

%Define the parameters that I'm not changing frequently:
N = 1;
%Don't need alpha since it's set to 5 in your lCDF function already

F(1) = gj.*(lCDF(r(1) + L/2) - lCDF(r(1) - L/2)) + ...
        (gj1 - gj).*(lCDF(2.*r(1)) - lCDF(0)) + ...
        (gj2 - gj1).*(lCDF(r(1) + r(2)) - lCDF(r(1) - r(2))) - j*N/m;
F(2) = gj.*(lCDF(r(2) + L/2) - lCDF(r(2) - L/2)) + ...
        (gj1 - gj).*(lCDF(r(2) + r(1)) - lCDF(r(2) - r(1))) + ...
        (gj2 - gj1).*(lCDF(2.*r(2)) - lCDF(0)) - (j+1)*N/m;
    
end