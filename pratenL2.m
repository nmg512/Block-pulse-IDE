function nL2 = pratenL2(L, p, gvals, r)

% Calculates n(L/2) for a given p-level equilibrium
% Uses equation 51 (re-indexed appropriately) evaluated at x = L/2

% Function takes arguments of:
% L: domain length;
% p: number of levels in the equilibrium;
% gvals: appropriate growth coefficients g_i (vector);
% r: the spatial-threshold fixed-points r_i (vector).

%First term in n(L/2)
term1 = gvals(1).*(lCDF(L) - lCDF(0));
full = term1;
%Calculate the next terms in the sum and add them successively
for j = 1:p-1
    nextterm = (gvals(j+1) - gvals(j)).*(lCDF(L/2 + r(j)) - lCDF(L/2 - r(j)));
    full = full + nextterm;
end
%Get n(L/2)
nL2 = full;

end