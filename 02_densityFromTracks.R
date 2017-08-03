

rm(list=ls())

setwd('C:/Users/Clara PERON/Documents/PELAGIC2/1-Methods/Shearwaters/gps_integration-master/markov_density/')

require(Matrix) # to create sparse matrices
library(matlab)
library(raster)
library(RColorBrewer)
library(CLmapping)
library(fields)
library(gdata)
source("markovdens.R")
source("secondtest.R")

# Load data
  load('Shearw_data.Rdata')

  # Format data
    names(data) <- c('x', 'y', 'id', 'ind')
    head(data)
    
    # Create a raster grid
      r <- raster(xmn=3.5, ymn=42.5, xmx=6.5, ymx=43.5, res=0.05)
      values(r) <- 1:(dim(r)[1]*dim(r)[2])
      
    # Sum the number of locs per grid cell
      datasp <- data
      datasp$occ <- 1
      coordinates(datasp) <- ~x+y
      rtracks <- rasterize(datasp, r, 'occ', fun=sum)
      mxv <- max(values(rtracks), na.rm=T) # determine the colony loc
      
      plot(rtracks)
    
    # Set the location of the colony to 0
      rtracks[rtracks!=mxv] <- 1
      rtracks[rtracks==mxv] <- 0
      plot(rtracks)

    # Extract the colony location info
      data$qqq <- raster::extract(rtracks, datasp)

    # Extract the cell number  
      data$cc <- raster::extract(r, datasp)

    # Convert raster r into grid  
      grid <- as.data.frame(r, xy=T)
      names(grid) <- c('xgrid', 'ygrid', 'cc')

    # merge the two datasets to add the NA where no locs in the grid
      datac <- merge(data, grid, by='cc', all.y=T)

      datac <- drop.levels(datac)
      
      datan <- datac[is.na(datac$ind)==F,]

      save(datac, file='Shearw_test.Rdata')

# -- END FORMATTNG ------------------------------------------------------------------------------------------------------

      # Calculate Track-based and Markov densities
        MD = markovdens(datac$cc, datac$qqq)
      
        ## ERROR : Error in eigen(Conj(t(uu))) : infinite or missing values in 'x' 

      # perform test for second-order dependence and print statistics
        test.out = secondtest(cc, da=colo, extracell=TRUE)
        print(test.out)

      # 1) Create a raster with the Track based densities
        tab <- data.frame(nlocs=values(rtracks))
        tab <- tab[is.na(tab$nlocs)==F,]
        tdens <- MD[[2]] 
        tdens <- tdens[tdens>0,]
        r1 <- r
      
        values(r1) <- MD[[2]]
        r1[r1==0] <- NA
        plot(r1, col=tim.colors(20))
      
      # 2) Create a raster with the Markov based densities  
        r2 <- raster(nrows=numdiv, ncols=numdiv, xmn=min(xm), xmx=max(xm), ymn=min(ym),  ymx=max(ym))
        values(r2) <- flipud(t(matrix(MD[[1]], numdiv, numdiv)))
        r2[r2==0] <- NA
        plot(r2, col=tim.colors(20))
      
        ## plot results
        cols <- colorRampPalette(c('#fef0d9','#fdcc8a','#fc8d59','#e34a33','#b30000'))
        zlim = c(0, max(MD[[1]], MD[[2]]))
        
        par(mfrow=c(1,2), mar=c(1,1,1,1), omi=c(0.1,0.2,0.1,0.1), pty='s')
        #points(dat$x[dat$id==1], dat$y[dat$id==1], type='l', lwd=0.1, col='red')
        #sapply(1:max(dat$id), function(i) lines(y~x, dat, subset=id==i, lwd=0.1, col='red'))
        plot(log(r1+1), col=cols(15), main="Track density", zlim=zlim, xlim=xm, ylim=ym, add=T)
        
        #points(dat$x[dat$id==1], dat$y[dat$id==1], type='l', lwd=0.1, col='red')
        #sapply(1:max(dat$id), function(i) lines(y~x, dat, subset=id==i, lwd=0.1, col='red'))
        plot(log(r2+1), col=cols(15), main="Markov densities", zlim=zlim, xlim=xm, ylim=ym, add=T)
      
