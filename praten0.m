function n0 = praten0(L, p, gvals, r)

% Calculates n(0) for a given p-level equilibrium
% Uses equation 51 (re-indexed appropriately) evaluated at x = 0

% Function takes arguments of:
% L: domain length;
% p: number of levels in the equilibrium;
% gvals: appropriate growth coefficients g_i (vector);
% r: the spatial-threshold fixed-points r_i (vector).

%First term in n(0)
term1 = gvals(1).*(lCDF(L/2) - lCDF(-L/2));
full = term1;
%Calculate the next terms in the sum and add them successively
for j = 1:p-1
    nextterm = (gvals(j+1) - gvals(j)).*(lCDF(r(j)) - lCDF(-r(j)));
    full = full + nextterm;
end
%Get n(0)
n0 = full;

end