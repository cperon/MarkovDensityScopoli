The example R code provided calculates track-based and Markov chain relative densities from simulated tracking data 
(see main text for details on the simulated data), performs a statistical test for second order dependence, 
and plots a comparison between the track-based and Markov chain relative densities.

The code requires the Matrix R package.
 

To evaluate the code, ensure that all 3 R source files and the data files are in an R working directory. 
Source in the mdensity.R function and, for example, type:

out = mdensity("rand_uniformR.txt")

Calculated relative densities 
are output to the "out" object, results of the test for second-order dependence are printed to the command window, and a 3-panel plot is generated.