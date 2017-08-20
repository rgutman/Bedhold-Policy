#loading the data that was estimated in parallel for hospitalization
allMat = read.csv("simSynth_1_50.csv",header=TRUE)
tempMat = read.csv("simSynth_51_100.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_101_150.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_151_200.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_201_250.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_251_300.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_301_350.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_351_400.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynth_401_450.csv",header=TRUE)
allMat = rbind(allMat, tempMat)

#treated facilities
allMatT = subset(allMat,allMat[,3]==1)

#state mean and variance per quarter
varStateEff = sqrt(apply(allMatT[,4:11],1,var)*7/8)
meanStateEff = apply(allMatT[,4:27],2,mean) 

#overall difference
hospAll = sum(allMatT[,12:27]) / (16*dim(allMatT)[1])

#control facilities
allMatC = subset(allMat,allMat[,3]==0)

#performing "placebo test" for overall hospitalization on a sample (to reduce complexity of computation)
placeboVal = NULL
placeboAllVal = NULL
for(numIter in 1:100)
{
  vecOrd = NULL
  for(i in 1:length(allMatT[,1]))
  {
    vecOrd = c(vecOrd,sample(1:40,40))
  }
  for(i in 1:40)
  {
   allmatSub = subset(allMatC, vecOrd == i)
   placeboVal = rbind(placeboVal, apply(allmatSub[,4:27], 2, mean))
   placeboAllVal = c(placeboAllVal,sum(allmatSub[,12:27]) / (16*dim(allmatSub)[1]))
  }
}

#placebo quantiles and median
med = intEst =  apply(placeboVal,2,median)
intEst =  apply(placeboVal,2,quantile,c(0.025,0.975))

#overall placebo test one sided p-value
hospPlacebo = quantile(placeboAllVal,c(0.025,0.975))
hospPval = length(which(placeboAllVal >= hospAll)) / length(placeboAllVal)
###end placebo test

#creating graph for hospitalization as in the paper for simulated data
pdf("Fig1SynHospSim.pdf")
plot(1:24, meanStateEff, type="l", lwd=4, ylim=c(-0.04,0.04), xaxt='n',xlab = "Year (quarter)", ylab="Changes in Hospitalization Rate (per patient month)")
axis(1, at=seq(1,24,4), tick = seq(1,24,4), labels=1999:2004)
abline(v=8)
abline(h=0)

lines(1:24, med, lty=1,lwd=4,col="gray")
lines(1:24, intEst[1,], lty=2,lwd=4,col="gray")
lines(1:24, intEst[2,], lty=2,lwd=4,col="gray")
dev.off()


#####################Mortality###############################################
#loading the data that was estimated in parallel
allMat = read.csv("simSynthD_1_50.csv",header=TRUE)
tempMat = read.csv("simSynthD_51_100.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_101_150.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_151_200.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_201_250.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_251_300.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_301_350.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_351_400.csv",header=TRUE)
allMat = rbind(allMat, tempMat)
tempMat = read.csv("simSynthD_401_450.csv",header=TRUE)
allMat = rbind(allMat, tempMat)

#treated facilities
allMatT = subset(allMat,allMat[,3]==1)
varStateEff = sqrt(apply(allMatT[,4:11],1,var)*7/8)

meanStateEff = apply(allMatT[,4:27],2,mean) 
mortAll = sum(allMatT[,12:27]) / (16*dim(allMatT)[1])

#control facilities
allMatC = subset(allMat,allMat[,3]==0)


#performing "placebo test" using samples from the overall values
placeboVal = NULL
placeboAllVal = NULL
for(numIter in 1:100)
{
  vecOrd = NULL
  for(i in 1:length(allMatT[,1]))
  {
   vecOrd = c(vecOrd,sample(1:40,40))
  }

  for(i in 1:40)
  {
   allmatSub = subset(allMatC, vecOrd == i)
   placeboVal = rbind(placeboVal, apply(allmatSub[,4:27], 2, mean))
   placeboAllVal = c(placeboAllVal,sum(allmatSub[,12:27]) / (16*dim(allmatSub)[1]))
  }
}

#placebo median and quantiles for each quarter
med = intEst =  apply(placeboVal,2,median)
intEst =  apply(placeboVal,2,quantile,c(0.025,0.975))

#overall calculation one sided  p-value in mortPval
mortPlacebo = quantile(placeboAllVal,c(0.025,0.975))
mortPval = length(which(placeboAllVal >= mortAll)) / length(placeboAllVal)

#generating graphs similar to the one in the paper for simulated data
pdf("Fig1SynMortSim.pdf")
plot(1:24, meanStateEff, type="l", lwd=4, ylim=c(-0.04,0.04), xaxt='n',xlab = "Year (quarter)", ylab="Changes in Mortality Rate (per patient month)")
axis(1, at=seq(1,24,4), tick = seq(1,24,4), labels=1999:2004)
abline(v=8)
abline(h=0)

lines(1:24, med, lty=1,lwd=4,col="gray")
lines(1:24, intEst[1,], lty=2,lwd=4,col="gray")
lines(1:24, intEst[2,], lty=2,lwd=4,col="gray")
dev.off()