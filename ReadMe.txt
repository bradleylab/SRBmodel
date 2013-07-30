These files support the sulfate-reduction model described in 

A.S. Bradley, W.D. Leavitt, D.T. Johnston, 2011. Revisiting the dissimilatory sulfate-reduction pathway
Geobiology  vol 9, page 446-457

For questions or comments on these files or the paper, please contact Alex Bradley: abradley@eps.wustl.edu

=============
These files are in the programming language R. R is a statistical computing language. It was developed for statisticians, but can perform many of the same functions of MATLAB. Unlike MATLAB however, R is open source. 

R can be downloaded from this site:
http://cran.r-project.org

There is now also an integrated development environment for R, which may be downloaded here:
http://www.rstudio.org/ 
=============

Using the files:

There are two files
1) SRBparameters.R
2) SRB070511.R

The first file contains all the input parameters for the model. The second file (SRB070511.R) contains the numerical calculations. In most instances, the second file doesn't need to be opened or modified. However, both files should be in the same directory.

First, change the working directory of your R program to the directory containing these two files.

The parameters in SRBparameters.R are outlined in the paper. 
There are 5 sets of parameters


#Fractionation factors
This first set of parameters allows you to set the fractionation factor for each of the 21 reactions in the model. These are listed as a0 - a18, plus reverse fractionation factors for a0-a2

#Reference values
The second set of parameters sets the #reference values. These are the lambda values that determine for which minor isotope system you wish to make calculations.  

#calculation parameters
The value "RecDepth" determines the recursion depth of the program. This is used for two purposes. i) it is the recursion depth for the calculation of the fluxes. ii) it is the power to which the defining matrix will be raised for the calculation of the isotopic composition of each parameter. We have generally found that 1e6 is a good number here, as it is high enough to achieve stable solutions, but low enough for computation times to be reasonable. 

#Plot parameters
These parameters simply set the bounds of the plot that is output by the program. 
plotspec: determines which species will be plotted. Enter 11 for sulfide. 2 for sulfate. Codes are given in the annotation. 
 
x.axis.min= x axis minimum delta value
x.axis.max = x axis maximum delta value
y.axis.min = y axis minimum DELTA value (e.g. CAP delta 33S or 36S
y.axis.max = y axis maximum DELTA value (e.g. CAP delta 33S or 36S

#run code
The final line of this file calls the code that runs the program. 
=============

When parameters are set as you wish, simply copy the contents of SRBparameters.R, paste into the R console, and press RETURN. The program will take a few seconds to run. Run times will depend on your computer. We ran these files on a 2.53 GHz MacBook Pro, and runtimes were generally 12-15 seconds. 