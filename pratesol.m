function nx = pratesol(L, p, gvals, r, x)

% Calculates n(x) for a given p-level equilibrium
% Uses equation 51 (re-indexed appropriately)

% Function takes arguments of:
% L: domain length;
% p: number of levels in the equilibrium;
% gvals: appropriate growth coefficients g_i (vector);
% r: the spatial-threshold fixed-points r_i (vector);
% x: grid of the spatial domain (vector).

%First term in n(x)
term1 = gvals(1).*(lCDF(x + L/2) - lCDF(x - L/2));
full = term1;
%Calculate the next terms in the sum and add them successively
for j = 1:p-1
    nextterm = (gvals(j+1) - gvals(j)).*(lCDF(x + r(j)) - lCDF(x - r(j)));
    full = full + nextterm;
end
%Get n(x)
nx = full;

end