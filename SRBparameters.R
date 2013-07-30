#####################
#SRB Network model. For reference, please see:
#Revisiting the dissimilatory sulfate-reduction pathway
#A.S. Bradley*, W.D. Leavitt, D.T. Johnston, 
#Geobiology, in revision
#copyright 2011
#*for questions/comments email to: bradley@fas.harvard.edu
#####################

#Fractionation factors
a0 = 1.003					
a1 = 1
a2 = .975
a3 = .947
a4 = 1
a5 = 1
a6 = 1
a7 = 1.000
a8 = 1
a9 = 1.000
a10 = 1
a11 = 1.000
a12 = 1
a13 = 1
a14 = 1
a15 = 1.000
a16 = 1.000
a17 = 1.000
a18 = 1.000
a0r = 1.000
a1r = 1.000
a2r = 1.000

#reference values
lambda3x = .5147		#lambda value for 3xS. Apply for 33S or 36S as necessary
lambda3xref = .5147		#lambda value for reference line


#fluxes
#rev0=0.0001			#Reversibility at j0/j0r. Equivalent to DJ07 f3. Range 0-0.99
rev1 = 0.99				#Reversibility at j1/j1r. Should be full value (0.99) to approximate DJ07 model. Range 0-0.99
#rev2=.99 				#Reversibility at j2/j2r. Approximately equivalent to DJ07 f5. Range 0-0.99
dsrs1 = 0.0				#dsr1 def j6. Range 0-1
dsrs2 = 0.0				#dsr2 def j8. Range 0-1
leakSO3=0.0  			#Cell leakiness for sulfite. Range 0-1
leakS3O6=0.0 			#Cell leakiness for trithionate. Range 0-1
leakS2O3=0.0			#Cell leakiness for thiosulfate. Range 0-1
TR1 = 0.0				#proportion of trithionate reduction through TR-1 (as opposed to TF)

#calculation parameters
RecDepth=1e6			#recursion depth.

#Plot parameters - species to plot and limits of graph
plotspec = 11			#sulfide = 11, sulfite = 4, APS = 3, SO4in = 2
						#S3O6-1=7, S3O6-2=8, S2O3-1=9, S2O3-2=10
x.axis.min= -30
x.axis.max= 0
y.axis.min= 0
y.axis.max= .15

#run code
source("SRB070511.r")		