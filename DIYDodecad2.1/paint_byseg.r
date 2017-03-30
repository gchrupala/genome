paint_byseg <- function(calc, chr, region = c(), src='byseg.txt', width=1366, height=768, lwd=2, choice=c(), top=-1,bottom=-1, tofile=F)
{
   X <- read.table(src,skip=1)
   K <- dim(X)[2]-5
   Elements <- c(as.matrix(read.table(paste(calc,".txt",sep=""))))
   
   medseg <- (X[,2]+X[,3])/2  
   
   if (length(region)==0) {
      Y <- X[X[,1]==chr, 6:(5+K)]     
      which_ones <- (X[,1]==chr)
      POS <- (X[which_ones,2] + X[which_ones,3])/(2000000)
   }
   else {
      region_start <- round(region[1]*1000000)
      region_end <- round(region[2]*1000000)
      which_ones <- (X[,1]==chr & medseg>=region_start & medseg <= region_end)
      Y <- X[which_ones, 6:(5+K)]     
      POS <- medseg[which_ones]/1000000
   }
   
   
      N <- dim(Y)[1]
      if (N<1) {
         return("No SNPs in range")
      }
      
      
      dummy <- vector(length=N-1)
      dummy[1:(N-1)] <- 0
      dummy <- c(max(Y), dummy, 0)
      toplot <- 1:K

      MEANPROPS <- vector(length=K)
      for (i in 1:K) {
        MEANPROPS[i] = mean(Y[,i])
      }
   

      if (length(choice)>0) {
        toplot <- c()
        for (i in 1:K) {
          if (sum(Elements[i]==choice)>0)
            toplot <- c(toplot,i)
        }
      }
      else if (top>0) {
        toplot = order(MEANPROPS, decreasing=T)[1:top]
      }
      else if (bottom>0) {
        toplot = order(MEANPROPS, decreasing=F)[1:bottom]
      }
      if (tofile==T) {
        png(filename=paste("Chromosome",chr,'.png',sep=''), width=width,height=height)
        par(mar=c(12,6,6,7))
      }
      plot(c(POS-(max(POS)-min(POS))/4,max(POS)), dummy , col="#ffffff", xlab="Position (Mb)", ylab="Proportion (%)", main=paste("Chromosome #",chr))
         legend("topleft",Elements[toplot], col=rainbow(length(toplot)), pch=20 )

      count = 0      
      for (j in toplot) {
        count = count+1
 
        lines(POS,Y[,j],col=rainbow(length(toplot))[count], lwd=lwd)
        points(POS,Y[,j],col=rainbow(length(toplot))[count])
      }
      if (tofile==T) {
        dev.off()
      }
      
}
