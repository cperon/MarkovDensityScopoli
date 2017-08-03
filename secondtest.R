secondtest = function(cc, da, extracell){

## Whitehead & Jonsen
##
## Test for second-order dependence. 
##
##	Created 20/09/2012
##
## Arguments: 
##		cc = vector giving the sequence of cells, numbered from 1, ...
##		da = vector (same length as cc) of ones, except zero when tracks start
##		extracell = is an extra cell used to represent outside the study area


numc = max(cc)

if(extracell){
	xxx = which((cc[1:(length(cc)-2)] < numc) & (cc[2:(length(cc)-1)] < numc) &
		 (cc[3:length(cc)] < numc) & da[2:(length(da)-1)] & da[3:length(da)])
		 }
else {
	xxx = which(da[2:(length(da)-1)] & da[3:length(da)])
	}		 


## triplets in main study area
returns = sum(cc[xxx] == cc[xxx+2])	# return to same start
totaltrips = length(xxx) # total triplets
tqq = sparseMatrix(i=cc[xxx], j=cc[xxx+1], x=1, 
	dims=rep(numc-extracell, 2)) # triplets starting off in i,j
tmat = sparseMatrix(i=cc[1:(length(cc)-1)], j=cc[2:length(cc)], x=da[2:length(da)],
	dims=rep(numc,2))
	
if(extracell){
	tmat = tmat[1:(dim(tmat)[1]-1), 1:(dim(tmat)[2]-1)]
	}
## Calculate statistics
uu = tmat / (0.000001 + rowSums(tmat) * rep(1, dim(tmat)[2]))
exreturns = sum(colSums(tqq*t(uu)))
chi2 = ((returns - exreturns)^2) * (exreturns^-1 + (totaltrips-exreturns)^-1)
G = 2 * (returns * log(returns/exreturns) + (totaltrips-returns) *
	log((totaltrips-returns)/(totaltrips-exreturns)))	
	
print('Results of test for 2nd order dependence')
sprintf("no. ABA's %5.1f; expected no. ABA's %5.1f; chi-squared =%5.2f with 1 df (P=%6.4f); G =%5.2f with 1 df (P=%6.4f)", returns, exreturns, chi2, 1-pchisq(chi2,1), G, 1-pchisq(G,1))

}	