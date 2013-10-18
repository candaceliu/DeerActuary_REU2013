require(lattice)
results = read.csv('allData.csv')

results$sumx  = as.double(results$sumx)
results$sumx2 = as.double(results$sumx2)
results$summ  = as.double(results$summ)
results$summ2 = as.double(results$summ2)

results$mx = results$sumx/results$N
results$sx = ((results$sumx2-(results$sumx^2)/results$N)/(results$N-1))
results$sx[results$sx<0] = 0.0
results$sx = sqrt(results$sx)

alphaVals = unique(results$alpha)
P     = unique(results$P)

maxM = max(results$mx)
minM = min(results$mx)
maxS = max(results$sx)
minS = min(results$sx)

filename = "fullParameters.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,pch=0,xlim=c(430000,529333),ylim=c(0.0,0.15),
         main="Parameter Space",xlab="P",ylab=expression(alpha))
    for(P in seq(430000,530000,by=(530000-430000)/150))
      {
        for (a in seq(0.0,0.15,by=0.15/50))
          {
            points(P,a,pch=18,cex=0.6)
          }
      }
    dev.off()
  }


filename = "meanDeerByAlpha.png"
if(!file.exists(filename))
  {
    # Display the mean deer populations for different values of P
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minM,maxM),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Mean of the Deer Population',
         main=expression(paste("Sample Mean of the Deer Population vs. ",alpha)))
    symbol = 0;
    legendSymbols = double(0)
    legendText = character(0)
    for (P in seq(430000,530000,by=15*(530000-430000)/150))
      {
        symbol = symbol+1
        legendSymbols = c(legendSymbols,symbol)
        legendText = c(legendText,paste("P=",P))
        points(results$alpha[results$P == P],results$mx[results$P == P],
               pch=symbol,col=symbol)
      }
    legend(0,20100,legendText,pch=legendSymbols,col=legendSymbols)
                                        #stop("a")
    dev.off()
  }


filename = "meanDeerByAlpha_boxplot.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minM,maxM),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Mean of the Deer Population',
         main=expression(paste("Sample Mean of the Deer Population vs. ",alpha)))
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        boxplot(results$mx[results$alpha == a],add=TRUE,at=a,boxwex=0.015)
      }
                                        #stop("a")
    dev.off()
  }


filename = "stdDevDeerbyAlpha.png"
if(!file.exists(filename))
  {
    # Display the std. dev of the deer populations for different values of P
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minS,maxS),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Standard Deviation of the Deer Population',
         main=expression(paste("Sample Standard Deviation of the Deer Population vs. ",alpha)))
    symbol = 0;
    legendSymbols = double(0)
    legendText = character(0)
    for (P in seq(430000,530000,by=(530000-430000)/150))
      {
        symbol = (symbol + 1)%%30
        legendSymbols = c(legendSymbols,symbol)
        legendText = c(legendText,paste("P=",P))
        points(results$alpha[results$P == P],results$sx[results$P == P],
               col=symbol,pch=symbol)
      }
                                        #stop("a")
    dev.off()
  }


filename = "stdDevDeerbyAlpha_boxplot.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minS,maxS),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Standard Deviation of the Deer Population',
         main=expression(paste("Sample Standard Deviation of the Deer Population vs. ",alpha)))
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        boxplot(results$sx[results$alpha == a],add=TRUE,at=a,boxwex=0.015)
      }
                                        #stop("a")
    dev.off()
  }

results$mF = results$summ/results$N*100.0;
results$sFx = ((results$summ2-(results$summ^2)/results$N)/(results$N-1))*10000.0;
results$sFx[results$sx<0.0] = 0.0
results$sFx = sqrt(results$sx)


maxF = max(results$mF)
minF = min(results$mF)
maxsF = max(results$sFx)
minsF = min(results$sFx)

levelplot(results$mF~results$alpha*results$P,contour=TRUE,
          xlab=expression(alpha),ylab='P',main="Sample Mean Fund Balance")
levelplot(results$sFx~results$alpha*results$P,contour=TRUE,
          xlab=expression(alpha),ylab='P',main="Sample Std. Dev. Fund Balance")

filename = "meanfundBalanceByP.png"
if(!file.exists(filename))
  {
    # Display the mean deer populations for different values of P
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minF,maxF),xlim=c(min(results$P),max(results$P)),col=0,pch=0,
         xlab="P",ylab='Sample Mean of the Fund Level',
         main=expression(paste("Sample Mean of the Deer Population vs. ","P")))
    symbol = 0
    legendSymbols = double(0)
    legendText = character(0)
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        symbol = (symbol + 1)%%30
        legendSymbols = c(legendSymbols,symbol)
        legendText = c(legendText,as.expression(bquote(alpha == .(a))))
        points(results$P[results$alpha == a],results$mF[results$alpha == a],
               col=symbol,pch=symbol)
      }
    legend(430000,5.6E7,legendText,pch=legendSymbols,col=legendSymbols)
                                        #stop("a")
    dev.off()
  }


filename = "meanFundBalanceByP_boxplot.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minF,maxF),xlim=c(min(results$P),max(results$P)),col=0,pch=0,
         xlab="P",ylab='Sample Mean of the Fund Level',
         main=expression(paste("Sample Mean of the Deer Population vs. ","P")))
    for (P in seq(430000,530000,by=10*(530000-430000)/150))
      {
        boxplot(results$mF[abs(results$P-P)<1.0],add=TRUE,at=P,boxwex=10000)
      }
                                        #stop("a")
    dev.off()
  }


filename = "meanFundBalanceByAlpha.png"
if(!file.exists(filename))
  {
    # Display the mean deer populations for different values of alpha
    png(filename,width=1024,height=1024)
    #cat("yo\n")
    plot(0,0,ylim=c(minF,maxF),xlim=c(min(results$alpha),max(results$alpha)+0.05),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Mean of the Fund Level',
         main=expression(paste("Sample Mean of the Deer Population vs. ",alpha)))
    symbol = 0
    legendSymbols = double(0)
    legendText = character(0)
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        symbol = (symbol + 1)%%30
        legendSymbols = c(legendSymbols,symbol)
        legendText = c(legendText,as.expression(bquote(alpha == .(a))))
        points(results$alpha[results$alpha == a],results$mF[results$alpha == a],
               col=symbol,pch=symbol)
      }
    par(xpd=TRUE)
    legend(0.16,5.6E7,legendText,pch=legendSymbols,col=legendSymbols)
    par(xpd=FALSE)
                                        #stop("a")
    dev.off()
  }


filename = "meanfundBalanceByAlpha_boxplot.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minF,maxF),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
         xlab=expression(alpha),ylab='Sample Mean of the Fund Level',
         main=expression(paste("Sample Mean of the Deer Population vs. ",alpha)))
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        boxplot(results$mF[results$alpha == a],add=TRUE,at=a,boxwex=0.01)
      }
                                        #stop("a")
    dev.off()
  }


filename = "sdFundBalanceByP.png"
if(!file.exists(filename))
  {
    # Display the standard dev. of the deer populations for different values of alpha
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minsF,maxsF),xlim=c(min(results$P),max(results$P)+11000),col=0,pch=0,
         xlab="P",ylab='Sample Standard Deviaion of the Fund Value',
         main=expression(paste("Standard Deviation of the Fund Value vs. ","P")),
         bty='L')
    symbol = 0
    legendSymbols = double(0)
    legendText = character(0)
    for (a in alphaVals[seq(1,length(alphaVals),by=3)])
      {
        symbol = (symbol + 1)%%30
        legendSymbols = c(legendSymbols,symbol)
        legendText = c(legendText,as.expression(bquote(alpha == .(a))))
        points(results$P[results$alpha == a],results$sFx[results$alpha == a],
               col=symbol,pch=symbol)
      }
    par(xpd=TRUE)
    legend(531000,50,legendText,pch=legendSymbols,col=legendSymbols)
    par(xpd=FALSE)
    dev.off()
  }



filename = "sdfundBalanceByP_boxplot.png"
if(!file.exists(filename))
  {
    png(filename,width=1024,height=1024)
    plot(0,0,ylim=c(minsF,maxsF),xlim=c(min(results$P),max(results$P)+11000),col=0,pch=0,
         xlab="P",ylab='Sample Standard Deviaion of the Fund Value',
         main=expression(paste("Standard Deviation of the Fund Value vs. ","P")),
         bty='L')
    for (P in seq(430000,530000,by=10*(530000-430000)/150))
      {
        boxplot(results$sFx[abs(results$P-P)<1.0],add=TRUE,at=P,boxwex=10000)
      }
    dev.off()
  }
