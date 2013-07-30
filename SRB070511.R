#####################
#SRB Network model. For reference, please see:
#A branching network for sulfate reduction: implications for sulfur isotope fractionation
#A.S. Bradley, D.T. Johnston, W.D. Leavitt
#Pub Name, Year, pages, etc....
#copyright 2010
#for questions/comments email to: bradley@fas.harvard.edu
#####################

#
#Variables are defined in file "SRBparameters.R" which calls this file. 

#END of variable section
#note: vector is 
#SO4out=d3xS[1] 
#SO4=d3xS[1][2] 
#APS=d3xS[3] 
#SO3=d3xS[4] 
#S=d3xS[5] 
#S2plus=d3xS[6] 
#S3O6S1=d3xS[7] 
#S3O6S2=d3xS[8] 
#S2O3S1=d3xS[9] 
#S2O3S2=d3xS[10] 
#H2S=d3xS[11] 
#H2Sout=d3xS[12]
##################################################################################################
#All changes in this file change the source code. All variables are defined in SRBparameters.
#####################
#FUNCTIONS
##
#function to convert an isotope ratio 34S/32S to a delta34S value
R2Delta <-function(R,CDT)
{
Delta=((R/CDT)-1)*1000	
return(Delta)
	}

#function to convert a delta34S value to an isotope ratio 34S/32S
Delta2R <-function(Delta,CDT)
{
Ratio=((Delta/1000)+1)*CDT	
return(Ratio)
	}
	
# function to calculate the nth power of the matrix
matPower <- function(X,n)
{
    if(n != round(n)) {
        n <- round(n)
        warning("rounding exponent `n' to", n)
    }
    if(n == 0)
        return(diag(nrow = nrow(X)))
    n <- n - 1
    phi <- X
    ## pot <- X # the first power of the matrix.
    while (n > 0)
    {
        if (n %% 2)
            phi <- phi %*% X
        if (n == 1) break
        n <- n %/% 2
        X <- X %*% X
    }
    return(phi)
}
		
##
#end FUNCTIONS
#####################

#####################
#Variable definitions
##

#CDT ratios. cdt is listed as its known value. cdt3x is an approximate value. 
#These are 'placeholder' numbers; they cancel out of the equations and the specific value set here is irrelevant. cdt is set at its measured value. cdt3x is an approximate value for 33S/32S
cdt = .045005			#CDT ratio for 34/32. 
cdt3x = 0.007379		#CDT ratio for 3x/32. 

k=0									#k = counter for solution field 
m=0
j0 = 1								#j0 is always defined as 1.

##initialize 3xS alphas
a0_3x = a0^lambda3x
a1_3x = a1^lambda3x
a2_3x = a2^lambda3x
a3_3x = a3^lambda3x
a4_3x = a4^lambda3x
a5_3x = a5^lambda3x
a6_3x = a6^lambda3x
a7_3x = a7^lambda3x
a8_3x = a8^lambda3x
a9_3x = a9^lambda3x
a10_3x = a10^lambda3x
a11_3x = a11^lambda3x
a12_3x = a12^lambda3x
a13_3x = a13^lambda3x
a14_3x = a14^lambda3x
a15_3x = a15^lambda3x
a16_3x = a16^lambda3x
a17_3x = a17^lambda3x
a18_3x = a18^lambda3x
a0r_3x = a0r^lambda3x
a1r_3x = a1r^lambda3x
a2r_3x = a2r^lambda3x

#Matrices to keep track of values at the rev0 and rev2 values
Mrev0=matrix(0,100,100)
Mrev2=matrix(0,100,100)
d34SH2S=matrix(0,100,100)
d3xSH2S=matrix(0,100,100)
##
#End of variable definition
#####################

#####################
#Begin main loop						       
##

#iterate across various amounts of branched flux
for(rev2 in seq(0,.99,.01)){
m=m+1
k=0
for(rev0 in seq(0,.99,.01)) {
	k=k+1
	#fluxes
	j0 = 1
	j0r=rev0*j0
	forwardflux=j0-j0r		
	j1=forwardflux/(1-rev1)
	j1r=rev1*j1
	j2=forwardflux/(1-rev2)
	j2r=rev2*j2
	j15=leakSO3*forwardflux		
	#if there is no side flux, or no leak of trithionate/thiosuflate then j3 is easy
	j3=forwardflux-j15							
	#otherwise solution is recursive; find explicit solution at recursion depth n
	n=RecDepth	
	if(((-1 - dsrs1*leakS2O3 - dsrs2*leakS2O3 + 
     dsrs1*dsrs2*leakS2O3 - 2*dsrs1*leakS3O6 + 
     dsrs1*leakS2O3*leakS3O6)*(-dsrs1*leakS2O3 - dsrs2*leakS2O3 + 
     dsrs1*dsrs2*leakS2O3 - 2*dsrs1*leakS3O6 + 
     dsrs1*leakS2O3*leakS3O6))!=0){
		j3=(dsrs1*leakS2O3 + dsrs2*leakS2O3 - dsrs1*dsrs2*leakS2O3 +    2*dsrs1*leakS3O6 - dsrs1*leakS2O3*leakS3O6 -    dsrs1*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n -    dsrs2*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n +    dsrs1*dsrs2*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n -    2*dsrs1*leakS3O6*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n +    dsrs1*leakS2O3*leakS3O6*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n -    dsrs1*leakS2O3*leakSO3 - dsrs2*leakS2O3*leakSO3 +    dsrs1*dsrs2*leakS2O3*leakSO3 - 2*dsrs1*leakS3O6*leakSO3 +    dsrs1*leakS2O3*leakS3O6*leakSO3 - (-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^    n*leakSO3 - dsrs1*leakS2O3*rev0 - dsrs2*leakS2O3*rev0 +    dsrs1*dsrs2*leakS2O3*rev0 - 2*dsrs1*leakS3O6*rev0 +    dsrs1*leakS2O3*leakS3O6*rev0 +    dsrs1*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n*rev0 +    dsrs2*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n*rev0 -    dsrs1*dsrs2*leakS2O3*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^    n*rev0 +    2*dsrs1*leakS3O6*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n*rev0 -    dsrs1*leakS2O3*leakS3O6*(-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^n*rev0 +    dsrs1*leakS2O3*leakSO3*rev0 + dsrs2*leakS2O3*leakSO3*rev0 -    dsrs1*dsrs2*leakS2O3*leakSO3*rev0 + 2*dsrs1*leakS3O6*leakSO3*rev0 -    dsrs1*leakS2O3*leakS3O6*leakSO3*rev0 + (-dsrs2*leakS2O3 +       dsrs1*(-2*leakS3O6 + leakS2O3*(-1 + dsrs2 + leakS3O6)))^    n*leakSO3*rev0)/((-1 - dsrs1*leakS2O3 - dsrs2*leakS2O3 +      dsrs1*dsrs2*leakS2O3 - 2*dsrs1*leakS3O6 +      dsrs1*leakS2O3*leakS3O6)*(-dsrs1*leakS2O3 - dsrs2*leakS2O3 +      dsrs1*dsrs2*leakS2O3 - 2*dsrs1*leakS3O6 +      dsrs1*leakS2O3*leakS3O6))}
		#end of explicit solution for j3. 
		#if statement preserves previous definition of j3 if no side flux or thionate leaks
		#See RecursiveFluxSolution.nb for derivation of the explicit equation for j3
	j6 = dsrs1*j3
	j7 = 2*j6
	j4 = j3-j6
	j8 = dsrs2*j4
	j5 = j4-j8
	j9 = j8
	j16 = j6*leakS3O6
	j10 = j6 - j16
	j11 = j10*(1-TR1)
	j17 = (j8+j10)*leakS2O3
	j12 = j10+j10*(1-TR1)
	j13 = j10+j8-j17
	j14 = j13
	j18 = j5 + j13

	
	
    ################# Ratio solver ######################
	#Define Matrices
	# matrix A is the matrix for linear difference equations in d34S
	# matrix A3x is the matrix for linear difference equations in d3xS
	
	A =matrix(0,12,12)	
	A3x =matrix(0,12,12)				# initialize each matrix to zeros
	
	A[1,1]=1	#SO4out = 1 (e.g. 1*cdt, where cdt is input vector)
	A3x[1,1]=1	#SO4out = 1 (e.g. 1*cdt, where cdt is input vector)
	
	A[2,1]=a0*j0/(a1*j1 + a0r*j0r)  
	A3x[2,1]=a0_3x*j0/(a1_3x*j1 + a0r_3x*j0r)  #SO4 dep on SO4out
	
	A[2,3] =a1r*j1r/(a1*j1 + a0r*j0r) 
	A3x[2,3] =a1r_3x*j1r/(a1_3x*j1 + a0r_3x*j0r) #SO4 dep on APS
	
	A[3,2]=a1*j1/(a2*j2 + a1r*j1r)	
	A3x[3,2]=a1_3x*j1/(a2_3x*j2 + a1r_3x*j1r)	#APS dep on SO4
	
	A[3,4]=a2r*j2r/(a2*j2 + a1r*j1r)	
	A3x[3,4]=a2r_3x*j2r/(a2_3x*j2 + a1r_3x*j1r)	#APS dep on SO3
	
	A[4,3]=a2*j2/(a2r*j2r + a3*j3 + a7*j7 + a9*j9 + a11*j11 + a15*j15) 
	A3x[4,3]=a2_3x*j2/(a2r_3x*j2r + a3_3x*j3 + a7_3x*j7 + a9_3x*j9 + a11_3x*j11 + a15_3x*j15) #SO3 dep on APS
	
	A[4,8]=a12*j12/(a2r*j2r + a3*j3 + a7*j7 + a9*j9 + a11*j11 + a15*j15) 
	A3x[4,8]=a12_3x*j12/(a2r_3x*j2r + a3_3x*j3 + a7_3x*j7 + a9_3x*j9 + a11_3x*j11 + a15_3x*j15) #SO3 dep on S3O6S2
	
	A[4,10]=a14*j14/(a2r*j2r + a3*j3 + a7*j7 + a9*j9 + a11*j11 + a15*j15) 
	A3x[4,10]=a14_3x*j14/(a2r_3x*j2r + a3_3x*j3 + a7_3x*j7 + a9_3x*j9 + a11_3x*j11 + a15_3x*j15) # SO3 dep on S2O3S2
	
	A[5,6]= a4*j4/(a5*j5 + a8*j8) 
	A3x[5,6]= a4_3x*j4/(a5_3x*j5 + a8_3x*j8) #S dep on Splus
	
	A[6,4]= a3*j3/(a6*j6 + a4*j4) 
	A3x[6,4]= a3_3x*j3/(a6_3x*j6 + a4_3x*j4) #S2plus dep on SO3
	
	if({(a10*j10+j16)!=0}){	
		A[7,6]= a6*j6/(a10*j10 + j16) 
		A3x[7,6]= a6_3x*j6/(a10_3x*j10 + j16)}#S3O6S1 dep on S2plus
	
	if({(a12*j12 + 2*j16)!=0}){	
		A[8,4] = a7*j7/(a12*j12 + 2*j16) 
		A3x[8,4] = a7_3x*j7/(a12_3x*j12 + 2*j16)} #S3O6S2 dep on SO3
	
	if({(a13*j13 + j17)!=0}){	
		A[9,5] = a8*j8/(a13*j13 + j17)	
		A3x[9,5] = a8_3x*j8/(a13_3x*j13 + j17) #S2O3S1 dep on S
		
		A[9,7]=a10*j10/(a13*j13 + j17)
		A3x[9,7]=a10_3x*j10/(a13_3x*j13 + j17)} #S2O3S1 dep on S3O6S1
	
	if({(a14*j14 + j17)!=0}){
		A[10,4]= (a9*j9+a11*j11)/(a14*j14 + j17) 
		A3x[10,4]= (a9_3x*j9+a11_3x*j11)/(a14_3x*j14 + j17)}#S2O3S2 dep on SO3
		
	if({(a14*j14 + j17)!=0}){
		A[10,8]= (TR1*j10)/(a14*j14 + j17) 
		A3x[10,8]= (TR1*j10)/(a14_3x*j14 + j17)}#S2O3S2 dep on S3O6S2
	
	A[11,5] = a5*j5/(a18*j18)  
	A3x[11,5] = a5_3x*j5/(a18_3x*j18)  #H2S dep on S
	
	A[11,9] = a13*j13/(a18*j18) 
	A3x[11,9] = a13_3x*j13/(a18_3x*j18) #H2S dep on S2O3 S1
	
	A[12,11] = a18  
	A3x[12,11] = a18_3x #H2S out is end member (i.e. boundary condition)

	#initialize starting values = 0 permil
	R34S= t(matrix(cdt,1,12))
	R3xS= t(matrix(cdt3x,1,12))	

	#find the nth power of the matrix for the RecDepth
	An=matPower(A,RecDepth)
	An3x=matPower(A3x,RecDepth)

	#Solve for output Ratios and assign to d34S and d3xS
	if(k==1){
		d34S = R2Delta(An%*%R34S,cdt)
		d3xS = R2Delta(An3x%*%R3xS,cdt3x)}else {
		d34S = cbind(d34S,R2Delta(An%*%R34S,cdt))
		d3xS = cbind(d3xS,R2Delta(An3x%*%R3xS,cdt3x))}	
    ################# end Ratio solver ######################

	Mrev0[m,k]=rev0
	Mrev2[m,k]=rev2
	
	}			# end of the k increment loop to solve across variable reversibility
		d34SH2S[m,]=d34S[plotspec,]	
		d3xSH2S[m,]=d3xS[plotspec,]
		
		}			# end m loop		

D3xSH2S=d3xSH2S - 1000*((1+d34SH2S/1000)^lambda3xref-1)

xlab.name=expression(paste(delta^ 34,"S"))
ylab.name=expression(paste(Delta^ 33,"S"))

plotlines=1
plot(d34SH2S[plotlines,],D3xSH2S[plotlines,],type='l', xlab=xlab.name, ylab=ylab.name, xlim=c(x.axis.min,x.axis.max),ylim=c(y.axis.min,y.axis.max))
plotlines=seq(10,90,9)
for (pcount in plotlines){
lines(d34SH2S[pcount,],D3xSH2S[pcount,])	
}
pcount=100
lines(d34SH2S[pcount,],D3xSH2S[pcount,])

plotlines=c(1,seq(10,100,10))
for (pcount in plotlines){
lines(d34SH2S[,pcount],D3xSH2S[,pcount])	
}



