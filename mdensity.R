mdensity = function(input){

##	Whitehead & Jonsen
##
##  Wrapper function to call markovdens function, perform test for second-order
##	 dependence, and plot a comparison of track-based density vs. Markov density. This
##	 code was created to use simulated data provided in Whitehead & Jonsen, but can 
##	 easily be modified to accomodate real tracking data
##
##	Created 20/09/2012
##
## Arguments: input = name of input data file, including file extension, in quotes
##
## Example:
##
##		source("mdensity.R")
##		out = mdensity("rand_uniformR.txt")

require(Matrix) # to create sparse matrices
source("markovdens.R")
source("secondtest.R")

## read data file
dat = read.table(input, header=FALSE)
names(dat) = c('x','y')
dat$id = rep(1:300, each=52)
dat = subset(dat, !is.na(x))
dat$ind = rep(1:51, 300)
dat = dat[,c(3,4,1,2)]

numdiv = 6
## analysis limits
xm = c(0, 1)
ym = c(0, 1)
qqq = dat$ind
qqq[qqq==1] = 0
qqq[qqq>1] = 1

## index location data to discrete grid cells
xc = ceiling(numdiv * ((dat[,3] - xm[1]) / (xm[2] - xm[1])))
yc = ceiling(numdiv * ((dat[,4] - ym[1]) / (ym[2] - ym[1])))
cc = xc + (yc-1) * numdiv
cc[which((xc>numdiv) | (xc<1) | (yc>numdiv) | (yc<1))] = numdiv^2+1

## calculate Track-based and Markov densities
out = markovdens(cc,qqq)

## perform test for second-order dependence and print statistics
test.out = secondtest(cc, da=qqq, extracell=TRUE)
print(test.out)

## plot results
cols = colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", 
	"brown"))
par(mfrow=c(1,3), mar=c(2,2,2,1), omi=c(0.1,0.2,0.1,0.1), pty='s')
plot(y~x, dat, subset=id==1, type='l', lwd=0.1, xlim=c(0,1), ylim=c(0,1), axes=FALSE, 
	xlab="",ylab="", main="Tracks")
box()
sapply(2:300, function(i) lines(y~x, dat, subset=id==i, lwd=0.1))
mtext(substr(input, 1, nchar(input)-5), 2, 2, cex=1)
# scale displayed values relative to maximum of track densities & markov densities
zlim = c(0, max(out[[1]],out[[2]]))
image(matrix(out$track.dens,6,6), col=cols(10), axes=FALSE, main="Track density", 
	zlim=zlim); box()
abline(v=seq(-0.1,1.1,l=7))
abline(h=seq(-0.1,1.1,l=7))
image(matrix(out$markov.dens,6,6), col=cols(10), axes=FALSE, main="Markov density", 
	zlim=zlim); box()
abline(v=seq(-0.1,1.1,l=7))
abline(h=seq(-0.1,1.1,l=7))

## return Markov density output
out
}