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

alpha = unique(results$alpha)
P     = unique(results$P)

maxM = max(results$mx)
minM = min(results$mx)
maxS = max(results$sx)
minS = min(results$sx)

plot(results$alpha[results$P == 430000],results$mx[results$P == 430000])


# Display the mean deer populations for different values of P
plot(0,0,ylim=c(minM,maxM),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
     xlab=expression(alpha),ylab='Sample Mean of the Deer Population',
     main=expression(paste("Sample Mean of the Deer Population vs. ",alpha)))
symbol = 0
for (plevel in P)
  {
    symbol = (symbol + 1)%%30
    symbol = 1
    points(results$alpha[results$P == plevel],results$mx[results$P == plevel],
           col=symbol,pch=symbol)
  }
#stop("yo!")

# Display the std. dev of the deer populations for different values of P
plot(0,0,ylim=c(minS,maxS),xlim=c(min(results$alpha),max(results$alpha)),col=0,pch=0,
     xlab=expression(alpha),ylab='Sample Standard Deviation of the Deer Population',
     main=expression(paste("Sample Standard Deviation of the Deer Population vs. ",alpha)))
symbol = 0
for (plevel in P)
  {
    symbol = (symbol + 1)%%30
    symbol = 1
    points(results$alpha[results$P == plevel],results$sx[results$P == plevel],
           col=symbol,pch=symbol)
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
          xlab=expression(alpha),ylab='P')
levelplot(results$sFx~results$alpha*results$P,contour=TRUE,
          xlab=expression(alpha),ylab='P')

# Display the mean deer populations for different values of P
plot(0,0,ylim=c(minF,maxF),xlim=c(min(results$P),max(results$P)),col=0,pch=0,
     xlab="P",ylab='Sample Mean of the Deer Population',
     main=expression(paste("Sample Mean of the Deer Population vs. ","P")))
symbol = 0
for (a in alpha)
  {
    symbol = (symbol + 1)%%30
    points(results$P[results$alpha == a],results$mF[results$alpha == a],
           col=symbol,pch=symbol)
  }


# Display the standard dev. of the deer populations for different values of alpha
plot(0,0,ylim=c(minsF,maxsF),xlim=c(min(results$P),max(results$P)),col=0,pch=0,
     xlab="P",ylab='Sample Mean of the Deer Population',
     main=expression(paste("Standard Deviation of the Deer Population vs. ","P")))
symbol = 0
for (a in alpha)
  {
    symbol = (symbol + 1)%%30
    points(results$P[results$alpha == a],results$sFx[results$alpha == a],
           col=symbol,pch=symbol)
  }
