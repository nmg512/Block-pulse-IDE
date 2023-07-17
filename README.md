# Block-pulse-IDE

Sample code to calculate equilibrium distributions for a ten-step block-pulse integrodifference equation (IDE) with an Allee growth function. This code sets up the block-pulse IDE and finds and plots all equilibria of the model.

A block-pulse IDE approximates the growth function of an IDE with a piecewise-constant function in order to create an analytically-tractable model. This method for analyzing IDEs using a block-pulse model is presented in the third chapter of the dissertation "Dynamics in Discrete Time: Successional Communities, Spatial Models, and Allee Effects" by N. M. Gilbertson, advised by M. Kot.

This code is intended as a sample only to provide a starting point for individuals wishing to explore block-pulse IDEs and their applications, and does not replicate any figures from that dissertation.

## Format

The associated code is in the form of Matlab .m files. The main file is `BP_equilibria.m`. This code builds the block-pulse IDE, sets up the different potential equilibria, and solves for and plots the valid equilibria of the block-pulse model.
The associated files `threeratesol.m`, `fourratesol.m`, etc. set up the systems of equations that are involved in solving for parameters used in the equilibrium expressions. The file `lCDF.m` contains the function for the Laplace cumulative distribution function. The files `pratesol.m`, `praten0.m`, and `pratenL2.m` provide general expressions for the full equilibrium distribution $n(x)$, maximum population density $n(0)$, and minimum population density $n(L/2)$ for an equilibrium with a given number of growth levels.

### License
This code is licensed under the terms of the MIT license. See the `LICENSE` file for details.
