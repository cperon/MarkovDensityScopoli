markovdens = function(cc, qqq, cellarea=NULL){

## Whitehead & Jonsen
##
## Markov density estimation from tracking data. This R function calculates 
##	and returns track-based and Markov relative densities.
##
##	Created 20/09/2012
##
## Arguments: 
##		cc = vector giving the sequence of cells, numbered from 1, ...
##		qqq = vector (same length as cc) of ones, except zero when tracks start
##		cellarea: vector of length max(cc) (i.e. number of cells) giving the area
##  		in each cell.  If the area of the final cell i given as NaN, then density
##  		is not estimated for this (exterior) cell	  

## number of cells
numcell = max(cc)
if(is.null(cellarea)){
	cellarea = c(rep(1, numcell-1), NA)
	}

extracell = as.numeric(is.na(cellarea[numcell]))

## create transition matrix
tmat = sparseMatrix(i=cc[1:length(cc)-1], j=cc[2:length(cc)], x=qqq[2:length(qqq)]
                    ,dims=c(numcell, numcell)   
                    )  
tmat = as.matrix(tmat) 

## remove unused cells
aa = which(colSums(tmat) > 0 | rowSums(tmat) > 0)  ## 20-03-2017 : clara replaced & by |
tmat = tmat[aa,aa]

## density estimation
uu = tmat/(rowSums(tmat) * rep(1, dim(tmat)[1]))
eig.out = eigen(Conj(t(uu))) # yields left eigenvalues and eigenvectors
eig.out$values = Re(eig.out$values) # take only the real part 
eig.out$vectors = Re(eig.out$vectors) # take only the real part

## Markov density
qj = which(eig.out$values == max(eig.out$values)) # index of eigenvalue 1
propta = eig.out$vectors[,qj] # eigenvecter assoicated w eigenvalue 1
propt = rep(0, numcell)
propt[aa] = propta

## correct for cellarea
propt = propt[1:(numcell - extracell)] / cellarea[1:(numcell - extracell)]

## standardize
propt = propt / sum(propt)

## track density estimate
orden = diag(sparseMatrix(i=cc,j=cc,x=1))
orden = orden[1:(numcell-extracell)] / cellarea[1:(numcell-extracell)]
orden = orden / sum(orden)

## return relative density estimates
list(markov.dens=propt, track.dens=orden) 
}