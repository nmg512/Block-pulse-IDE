function cdf_laplace = lCDF(x)

%Define parameters
alpha = 5;

%Output of function
cdf_laplace = (1/2).*(1 + sign(x).*(1 - exp(-alpha.*abs(x))));

end