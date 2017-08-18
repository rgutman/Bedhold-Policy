## Estimation of bed-hold policy
## approximate method using facility-level data
## Code for generating graphs and results with simulated data
## Roee Gutman 8/18/2017

require(mnormt)
require(lme4)
require(arm)
require(stringi)
require(reshape2)
require(MHadaptive)

# load R packages
library(mnormt)
library(lme4)
library(arm)
library(stringi)
library(reshape2)
library(MHadaptive)

#name of file that would be used for prior defenition 
fileLimits  = "fileLimitsSim30.csv"
 
#loading the simulated data
dataUse = read.csv("simulatedDataSetForPaper.csv",header=TRUE,row.names=1)


# Split states by bedhold policy
toBedHold = subset(dataUse,dataUse$BedHoldChange==1)
noBedHold = subset(dataUse,dataUse$BedHoldChange==0)

# noBedHoldF
#switching to rate log scale for no bed hold facilities
noBedHoldF = noBedHold
noBedHoldF[,1] = log(noBedHold[,1]/noBedHold[,3] + 0.0001) #acute hospitalization
noBedHoldF[,2] = log(noBedHold[,2]/noBedHold[,3] + 0.0001) #death
noBedHoldF[,3] = log(noBedHold[,3])                         #person months

facIDPos = unique(noBedHoldF$facid)
facIDPMatch = facIDPos
nNumPossMatch = length(facIDPMatch)
equalTo1 = rep(0,nNumPossMatch)
nNumBefore = 8

#creating facilities min and max hospitalization rate, mortality rate, and total person months file that would be used for prior defenition
matchLimitsMatrix = matrix(0, nNumPossMatch, 6)
for(j in 1:nNumPossMatch)
{
        XX = noBedHoldF[noBedHoldF$facid==facIDPMatch[j],1:3]
        matchLimitsMatrix[j,] = c(apply(XX,2,range))
}

write.table(matchLimitsMatrix, file=fileLimits,col.names = FALSE, row.names = FALSE,sep=",",qmethod="double")

#loading the file above which dictates the min/max on future values. This would be used for the prior that was used in the paper
limitTab = read.csv(fileLimits,header=FALSE)

#a function for not allowing param = (\beta_0, \beta_1) of regression to predict values
#beyond dMaxValue given mix and max value of x (minMax = (min value,max value))this implements the prior discussed in the paper
checkVariablesMinMax = function(param,minMax,dMaxValue)
{
   dMinVal = 0
   dMaxVal = 0
 	if (minMax[1] < 0 & minMax[2] < 0)
	{
		dMinVal = max((-param[1] + (dMaxValue)) / minMax[1], (-param[1] + (dMaxValue)) / minMax[2]);
		dMaxVal = Inf;
	}
	
	if (minMax[1] > 0 & minMax[2] > 0)
	{
		dMaxVal = min((-param[1] + (dMaxValue)) / minMax[1], (-param[1] + (dMaxValue)) / minMax[2]);
		dMinVal = -Inf;
	}

	if (minMax[2] > 0 & minMax[1] < 0)
	{
		dMinVal = (-param[1] + (dMaxValue)) / minMax[1];
		dMaxVal = (-param[1] + (dMaxValue)) / minMax[2];
	}

	if (param[2] > dMaxVal | param[2] < dMinVal)
	{
   return(Inf) #value is not allowed
  }
  return (0)  #value is allowed
}

#simple linear regression likelihood of yy on xx function with prior distribution that doesn't allow for 
#predicted values to be outside the range of min max
#parV = (\beta_0,\beta_1) sigma is the variance of residuals minMax=(max possible value of x, min possible value of x)
#maxVal is defined in advance as described in the paper
normLike = function(parV, yy, xx, sigma, minMax, maxVal)
{
    if(checkVariablesMinMax(parV,minMax,maxVal) ==0)
    {     
      beta0 <- parV[1]
      beta1 <- parV[2]
      # the mean of the regression
      meanV <- beta0 + beta1 * xx 
      # the log-likelihood of the whole dataset
      LL <- sum(dnorm(yy, mean=meanV, sd = sqrt(sigma), log = TRUE))
      return(LL)
    }
    return(-Inf)
}
#a function for sampling from the approximate posterior distribution, sigma can be sampled
#directly from an inverse chi square distribution, because of the priors beta is sampled using Metropolis Hastings algorithm
sampleCoeff = function(yy, xx, nNumSamp, minMax, dMaxValue) 
{
    numIter = 50
    numBurn = 10
    lmObj = lm(yy~xx)
    QR<-lmObj$qr
    df.residual<-lmObj$df.residual
    R<-qr.R(QR) ## R component
    coef2<-lmObj$coef
    Vb<-chol2inv(R) ## variance(unscaled)
    s2<-(t(lmObj$residuals)%*%lmObj$residuals)
    s2<-s2[1,1]/df.residual
 
    ## now to sample residual variance
    sigma<-df.residual*s2/rchisq(nNumSamp,df.residual)
    coefMat = NULL
    for(ii in 1:nNumSamp)#given sampled variance sampled regression coefficents using MH algorithm
    {
      mcmc_r = Metro_Hastings(normLike, coefficients(lmObj), prop_sigma = vcov(lmObj), iterations = numIter, burn_in = numBurn,  quiet =TRUE, yy=yy,xx=xx,minMax=minMax,maxVal=dMaxValue,sigma=sigma[ii])

      coefMat = rbind(coefMat,mcmc_r$trace[numIter - numBurn,])     
    }
   	return(list(coef=coefMat,sigma=sqrt(sigma)))
}
#likelihood function for simple Poisson regression of yy on xx with limiting prior distribution as described in the paper
#parV = (\beta_0,\beta_1) sigma is the variance of residuals minMax=(max possible value of x, min possible value of x)
#maxVal is defined in advance as described in the paper
poisLike = function(parV, yy, xx, offsetV, minMax, maxVal)
{
    if(checkVariablesMinMax(parV,minMax,maxVal) ==0)
    {     
      beta0 <- parV[1]
      beta1 <- parV[2]
  
      # the rate
      lambda <- exp(beta0 + beta1 * xx + offsetV)
      # the log-likelihood of the whole dataset
      LL <- sum(dpois(yy, lambda, log = TRUE))
      return(LL)
    }
    return(-Inf)
}

#function for sampling from the posterior distribution of hospitalization and mortality rates using MH algorithm
# and the function for the constrained posson likelihood described above
sampleCoeff2 = function(yy, xx, offsetV, nNumSamp, minMax, dMaxValue) 
{
  
  lmMod = glm(yy~xx+offset(offsetV),family=poisson(log))
  mcmc_r = Metro_Hastings(poisLike, coefficients(lmMod), prop_sigma = vcov(lmMod), iterations = 6*nNumSamp, burn_in = nNumSamp,  quiet = TRUE, yy=yy,xx=xx,offsetV=offsetV,minMax=minMax,maxVal=dMaxValue)

  mcmc_r<-mcmc_thin(mcmc_r,5)

  return(list(coef=mcmc_r$trace[1:nNumSamp,]))
}

######################This is just to check the function not required#########################
#sampB2 = sampleCoeff2(yValue2_b4, xValue2_obs[1:numPeriod], (log(yValue1_b4)), numImp,minMax=limitTab[matchfac[j],1:2],dMaxValue=log(6))
    
# use maximum likelihood to perform matching
#we have three models two poissons (mortality and hospitalization) and one log-normal for resident months)
#this is used as an approximation instead of the Bayesian matching described in the paper
calclik_3Pois = function(x, y, j) {
    
    #facility that did not change the policy
    xValue1 = x[,3]
    xValue2 = x[,1]
    xValue3 = x[,2]
    
    #facility that changed the policy
    yValue1 = y[,3]
    yValue2 = y[,1]
    yValue3 = y[,2]
    
    #resident months equation
    lmMod1 = lm(I(log(yValue1)) ~ xValue1)
    
    #acute hospitilization equation     
    lmMod2 = glm(yValue2 ~ xValue2 + offset(log(yValue1)), family = poisson())
    
    #mortality equation
    lmMod3 = glm(yValue3 ~ xValue3 + offset(log(yValue1)), family = poisson())

    #asymptotic variance/covariance of the parameters    
  mod1Mat = vcov(lmMod1)
  mod2Mat = vcov(lmMod2)
  mod3Mat = vcov(lmMod3)
  logLikeAll = 0
  #this is to take care of cases where one facility has almost constant rates or person months. Simply making sure that \beta_1 is not 0 - if it is the match would be invalid  
  if( sum(dim(mod3Mat) == c(2,2)) == 2 && sum(dim(mod2Mat) == c(2,2)) == 2 && sum(dim(mod1Mat) == c(2,2))==2)
  {
#if future predicted results are beyond plausible limits they are being given exteremly high distance - this is instead the implementation of the prior, so that not to encounter unreasonable values due to extrapolation
     if(checkVariablesMinMax(coef(lmMod1),minMax=limitTab[j,5:6],dMaxValue=log(2000))==Inf)
     {
      return (-Inf)
     }
     logLik1 = logLik(lmMod1)
     
     if(checkVariablesMinMax(coef(lmMod2),minMax=limitTab[j,1:2],dMaxValue=log(6))==Inf)
     {
      return (-Inf)
     }
     logLik2 = logLik(lmMod2)
    
      if(checkVariablesMinMax(coef(lmMod3),minMax=limitTab[j,3:4],dMaxValue=log(6))==Inf)
     {
      return (-Inf)
     }
    logLik3 = logLik(lmMod3)

     logLikeAll = logLik1 + logLik2 + logLik3
 } else
 {
  logLikeAll = -Inf;
 }
    return(logLikeAll)
}

#creating the matches
n.tobed = length(unique(toBedHold$facid))
n.nobed = length(unique(noBedHoldF$facid))
# calculate the loglikelihoods over all possible pairs (relatively slow process, but only need to be done once - this is why it was not optimized)
rSqMatchtoBedF_lik2 = matrix(0, nrow = n.tobed, ncol = n.nobed)
for (i in 1:n.tobed) {
#only choosing the first time point for each facility which are before the policy
    tobed.y = toBedHold[(1+24*(i-1)):(8+24*(i-1)),1:3]
    for (j in 1:n.nobed) {
        nobed.x = noBedHoldF[(1+24*(j-1)):(8+24*(j-1)),1:3]
        rSqMatchtoBedF_lik2[i,j] = calclik_3Pois(nobed.x, tobed.y,j)
    }
}

#saving it to a file because this takes very long so running it once
save(rSqMatchtoBedF_lik2, file="BedHold_DistMat_ApproxSim.Rdata")

# load distance matrix
#load("BedHold_DistMat_ApproxSim.Rdata")

#selecting the best match using a greedy algorithm for computational efficiency 
tempTab = rSqMatchtoBedF_lik2
matchfac = rep(-1,n.tobed)
numMatched = 0
for (i in 1:n.tobed) 
{
 value = which(tempTab == max(tempTab), arr.ind = TRUE)
 if(matchfac[value[1]] == -1)
 {
  matchfac[value[1]] = value[2]
  numMatched = numMatched +1
  print(numMatched)
  tempTab[,value[2]] = -Inf
  tempTab[value[1],] = -Inf
 }
 
}

matchedTo = NULL
n.fac = length(matchfac)
for(j in 1:n.fac) {
	k = matchfac[j]
    matchto = subset(noBedHoldF, noBedHoldF[,4] == unique(noBedHoldF[,4])[k])
    matchedTo = rbind(matchedTo, matchto)
}
##End greedy matching algorithm#################################

#partitioning observation in the state that changed policy to time periods before the change and time periodd after
toBedHoldNumBeforeF = subset(toBedHold, toBedHold[,5] == 0) 
toBedHoldNumAfterF = subset(toBedHold, toBedHold[,5] == 1)

xValue1 = matchedTo[,3] #resident Months facilities with no change of policy
xValue2 = matchedTo[,1] #acute hospitilization facilities with no change of policy
xValue3 = matchedTo[,2] #mortality facilities with no change of policy

yValue1 = toBedHoldNumBeforeF[,3] #resident Months facilities before change of policy
yValue2 = toBedHoldNumBeforeF[,1]  #acute hospitilization facilities before change of policy
yValue3 = toBedHoldNumBeforeF[,2]  #mortality facilities before change of policy

obsValue1 = toBedHoldNumAfterF[,3]  #resident Months facilities with after change of policy
obsValue2 = toBedHoldNumAfterF[,1]  #acute hospitilization facilities before change of policy
obsValue3 = toBedHoldNumAfterF[,2]  #mortality facilities after change of policy

set.seed(1234)    #for reproducibility reasons
# Imputation Arrays
numImp = 100
imputeArr1 = imputeArr2 = array(0, c(n.fac, length(obsValue1)/n.fac, numImp))
imputeArrPre1 = imputeArrPre2 = array(0, c(n.fac, length(obsValue1)/n.fac, numImp))
imputeArrVal1 = imputeArrVal2 = imputeArrVal3 = array(0, c(n.fac, length(obsValue1)/n.fac, numImp))

imputeArrValB1 = imputeArrValB2 = imputeArrValB3 = array(0, c(n.fac, length(yValue1)/n.fac, numImp))
#creating the imputations
for(j in 1:n.fac) {

#observed no policy for 24 periods    
    sel1 = (1+24*(j-1)):(24*j)
    xValue1_obs = xValue1[sel1]
    xValue2_obs = xValue2[sel1]
    xValue3_obs = xValue3[sel1]

#before any policy change outcomes for facilities that experienced the change    
    sel2 = (1+8*(j-1)):(8*j)
    yValue1_b4 = yValue1[sel2]
    yValue2_b4 = yValue2[sel2]
    yValue3_b4 = yValue3[sel2]

#after any policy change for facilities that experienced the change    
    sel3 = (1+16*(j-1)):(16*j)
    obsValue1_aftr = obsValue1[sel3]
    obsValue2_aftr = obsValue2[sel3]
    obsValue3_aftr = obsValue3[sel3]
    
    #number of period before and total
    numPeriod = length(yValue1_b4)
    totPeri = length(xValue1_obs)
    
    # 2 poisson regression models   and one linear regrssion for matched facilities
    lmMod1 = lm(I(log(yValue1_b4)) ~ xValue1_obs[1:numPeriod])
    lmMod2 = glm(yValue2_b4 ~ xValue2_obs[1:numPeriod] + offset(log(yValue1_b4)), family = poisson())
    lmMod3 = glm(yValue3_b4 ~ xValue3_obs[1:numPeriod] + offset(log(yValue1_b4)), family = poisson())
 
    
    # sample from posterior distribution
    sampB1 = sampleCoeff(log(yValue1_b4), xValue1_obs[1:numPeriod], numImp,minMax=limitTab[matchfac[j],5:6],dMaxValue=log(2000))
    sampB2 = sampleCoeff2(yValue2_b4, xValue2_obs[1:numPeriod], (log(yValue1_b4)), numImp,minMax=limitTab[matchfac[j],1:2],dMaxValue=log(6))
    sampB3 = sampleCoeff2(yValue3_b4, xValue3_obs[1:numPeriod], (log(yValue1_b4)), numImp, minMax=limitTab[matchfac[j],3:4],dMaxValue=log(6))
    
    for(i in 1:numImp){
 #imputing rates and resident months       
        lambda1 = rnorm(totPeri,mean=cbind(rep(1,totPeri),xValue1_obs) %*% sampB1$coef[i,],sd=sampB1$sigma[i])
        predValue1 = exp(lambda1)

        lambda2 = exp(cbind(rep(1,totPeri),xValue2_obs) %*% sampB2$coef[i,] + log(predValue1))
        predValue2 = rpois(totPeri, lambda2)
 
        lambda3 = exp(cbind(rep(1,totPeri),xValue3_obs) %*% sampB3$coef[i,] + log(predValue1))
 
        predValue3 = rpois(totPeri, lambda3)
        
#after policy change imputing missing potential control for state that changed policy        
        imputeArrVal1[j,,i] = predValue1[(numPeriod+1):totPeri]
        imputeArrVal2[j,,i] = predValue2[(numPeriod+1):totPeri]
        imputeArrVal3[j,,i] = predValue3[(numPeriod+1):totPeri]
        
        imputeArrPre1[j,,i] = predValue2[(numPeriod+1):totPeri] / predValue1[(numPeriod+1):totPeri]
        imputeArrPre2[j,,i] = predValue3[(numPeriod+1):totPeri] / predValue1[(numPeriod+1):totPeri]
     #calculating the difference for each facility for each post-policy time period   
        imputeArr1[j,,i] = obsValue2_aftr / obsValue1_aftr -  predValue2[(numPeriod+1):totPeri] / predValue1[(numPeriod+1):totPeri]
        imputeArr2[j,,i] = obsValue3_aftr / obsValue1_aftr -  predValue3[(numPeriod+1):totPeri] / predValue1[(numPeriod+1):totPeri]
#imputing before values as well (will be used in generating graphs)        
        imputeArrValB1[j,,i] = predValue1[1:numPeriod]
        imputeArrValB2[j,,i] = predValue2[1:numPeriod]/predValue1[1:numPeriod]
        imputeArrValB3[j,,i] = predValue3[1:numPeriod]/predValue1[1:numPeriod]
    }  # end of i
    print(j)
}  # end of j


allImpMB = NULL
allImpVB = NULL
allImpMBM = NULL
allImpVBM = NULL
allImpMBR = NULL
allImpVBR = NULL
#summarizing the imputation for all time before change happenes
#would be used to create graphs
for(i in 1:numImp){
#hospitilization
 mImpB = apply(imputeArrValB2[,,i],2,mean)
 vImpB = apply(imputeArrValB2[,,i],2,var) / length(imputeArrValB2[,1,i])
 
 allImpMB = rbind(allImpMB,mImpB)
 allImpVB = rbind(allImpVB,vImpB)
#mortality
 mImpBM = apply(imputeArrValB3[,,i],2,mean)
 vImpBM = apply(imputeArrValB3[,,i],2,var) / length(imputeArrValB3[,1,i])
 
 allImpMBM = rbind(allImpMBM,mImpBM)
 allImpVBM = rbind(allImpVBM,vImpBM)
#resident months
 mImpBR = apply(imputeArrValB1[,,i],2,mean)
 vImpBR = apply(imputeArrValB1[,,i],2,var) / length(imputeArrValB1[,1,i])
 
 allImpMBR = rbind(allImpMBR,mImpBR)
 allImpVBR = rbind(allImpVBR,vImpBR)
}

# function to calculate overall mean effect and variance within imputed dataset
#relies on multilevel model to adjust for year to year variability
summ.Arr = function(dat){
	# rearrange data 
	data = cbind(rep(1,dim(dat)[1]), dat)
	time = 1:dim(dat)[2]
	colnames(data) = c("subject", paste("t", time, sep = ""))
	data1 = melt(data, id.vars = "subject")
	data2 = data1[-c(1:dim(dat)[1]), ]
	colnames(data2)[1:2] = c("subject", "time")
	data3 = data2[ order(data2$subject), ]
	
	# lmer
	fit = lmer(value ~ 1 + (1|subject), data = data3)
	mean.fit  = as.numeric(fixef(fit))               # fixed effect coefficient
	coef.sigma.fit  = as.numeric(se.fixef(fit))      # standard error of fixed effect coefficient
	coef.sigma2.fit = as.numeric(se.fixef(fit)^2)    # variance of fixed effect coefficient
	
	return(list(mean = mean.fit, var = coef.sigma2.fit, sigma = coef.sigma.fit))
}


# multiple imputation
imputeMat1 = imputeMat2 = matrix(0, numImp, 3)
imputeMat1s = imputeMat2s = matrix(0, numImp, 3)
imputeMat1l = imputeMat2l = matrix(0, numImp, 3)

#for each facility finding the quantile that it was located at. 
#the quantile is similar for all quarters/years so the minimum was chosen
quantileNH = tapply(toBedHold$totbed,toBedHold$facid,min)

#summarizing total results over the entire post-policy period
for(i in 1:numImp){
	# rehospitalization rate
	summ1 = summ.Arr(imputeArr1[,,i])
	imputeMat1[i,1] = summ1$mean
	imputeMat1[i,2] = summ1$var
	imputeMat1[i,3] = summ1$sigma
	
	# mortality
	summ2 = summ.Arr(imputeArr2[,,i])
	imputeMat2[i,1] = summ2$mean
	imputeMat2[i,2] = summ2$var
	imputeMat2[i,3] = summ2$sigma
 
#Small NH 
  # rehospitalization rate
	summ1 = summ.Arr(imputeArr1[quantileNH == 1,,i])
	imputeMat1s[i,1] = summ1$mean
	imputeMat1s[i,2] = summ1$var
	imputeMat1s[i,3] = summ1$sigma
	
	# mortality
	summ2 = summ.Arr(imputeArr2[quantileNH == 1,,i])
	imputeMat2s[i,1] = summ2$mean
	imputeMat2s[i,2] = summ2$var
	imputeMat2s[i,3] = summ2$sigma

#large NH 
  # rehospitalization rate
	summ1 = summ.Arr(imputeArr1[quantileNH == 5,,i])
	imputeMat1l[i,1] = summ1$mean
	imputeMat1l[i,2] = summ1$var
	imputeMat1l[i,3] = summ1$sigma
	
	# mortality
	summ2 = summ.Arr(imputeArr2[quantileNH == 5,,i])
	imputeMat2l[i,1] = summ2$mean
	imputeMat2l[i,2] = summ2$var
	imputeMat2l[i,3] = summ2$sigma	
}

estAve1 = estAve2 = NULL    # overall estimate
estSd1  = estSd2  = NULL    # total standard error
estDef1 = estDef2 = NULL    # degrees of freedom
#results for overall
estAve1 = mean(imputeMat1[,1])
estSd1  = sqrt( mean(imputeMat1[,2]) + (1 + 1/numImp)*var(imputeMat1[,1]) )
estDef1 = (numImp - 1) * (1 + (numImp * mean(imputeMat1[,2])) / ((numImp + 1) * var(imputeMat1[,1])))^2

estAve2 = mean(imputeMat2[,1])
estSd2  = sqrt(mean(imputeMat2[,2]) + (1 + 1/numImp)*var(imputeMat2[,1]))
estDef2 = (numImp - 1) * (1 + (numImp * mean(imputeMat2[,2])) / ((numImp + 1) * var(imputeMat2[,1])))^2


upper1 = estAve1 + qt(0.975,estDef1)*estSd1
lower1 = estAve1 + qt(0.025,estDef1)*estSd1

upper2 = estAve2 + qt(0.975,estDef2)*estSd2
lower2 = estAve2 + qt(0.025,estDef2)*estSd2


estAve1s = estAve2s = NULL    # overall estimate
estSd1s  = estSd2s  = NULL    # total standard error
estDef1s = estDef2s = NULL    # degrees of freedom

#results for overall small NH
estAve1s = mean(imputeMat1s[,1])
estSd1s  = sqrt( mean(imputeMat1s[,2]) + (1 + 1/numImp)*var(imputeMat1s[,1]) )
estDef1s = (numImp - 1) * (1 + (numImp * mean(imputeMat1s[,2])) / ((numImp + 1) * var(imputeMat1s[,1])))^2

estAve2s = mean(imputeMat2s[,1])
estSd2s  = sqrt(mean(imputeMat2s[,2]) + (1 + 1/numImp)*var(imputeMat2s[,1]))
estDef2s = (numImp - 1) * (1 + (numImp * mean(imputeMat2s[,2])) / ((numImp + 1) * var(imputeMat2s[,1])))^2


upper1s = estAve1s + qt(0.975,estDef1s)*estSd1s
lower1s = estAve1s + qt(0.025,estDef1s)*estSd1s

upper2s = estAve2s + qt(0.975,estDef2s)*estSd2s
lower2s = estAve2s + qt(0.025,estDef2s)*estSd2s

estAve1l = estAve2l = NULL    # overall estimate
estSd1l  = estSd2l  = NULL    # total standard error
estDef1l = estDef2l = NULL    # degrees of freedom

#results for overall large NH
estAve1l = mean(imputeMat1l[,1])
estSd1l  = sqrt( mean(imputeMat1l[,2]) + (1 + 1/numImp)*var(imputeMat1l[,1]) )
estDef1l = (numImp - 1) * (1 + (numImp * mean(imputeMat1l[,2])) / ((numImp + 1) * var(imputeMat1l[,1])))^2

estAve2l = mean(imputeMat2l[,1])
estSd2l  = sqrt(mean(imputeMat2l[,2]) + (1 + 1/numImp)*var(imputeMat2l[,1]))
estDef2l = (numImp - 1) * (1 + (numImp * mean(imputeMat2l[,2])) / ((numImp + 1) * var(imputeMat2l[,1])))^2


upper1l = estAve1l + qt(0.975,estDef1l)*estSd1l
lower1l = estAve1l + qt(0.025,estDef1l)*estSd1l

upper2l = estAve2l + qt(0.975,estDef2l)*estSd2l
lower2l = estAve2l + qt(0.025,estDef2l)*estSd2l

#creating one table for overall results. This table is
allRes = rbind(c(lower1, upper1, lower2,upper2), c(lower1s, upper1s, lower2s,upper2s), c(lower1l, upper1l, lower2l,upper2l))

colnames(allRes) = c("lower 95 Hospitaliation","upper 95 Hospitaliation","lower 95 Mortality","upper 95 Mortality")
rownames(allRes) = c("All Facilities", "Small Facilities", "Large Facilities")

write.csv(round(allRes,4),"allPeriodRes.csv")
##############################################Creating Graphs########################################
#Hospitalization per month
checkVal = toBedHold[,1]/(toBedHold[,3])
#Mortality per month
checkValM = toBedHold[,2]/(toBedHold[,3])
#Residents per month
checkValR = toBedHold[,3]

#creating hospitalization, mortality and residents per state for year and quarter for change with policy change
newValChange = rep(1:24, length(unique(toBedHold$facid)))
forStateChange = tapply(checkVal, newValChange, mean)
forStateChangeM = tapply(checkValM, newValChange, mean)
forStateChangeR = tapply(checkValR, newValChange, mean)

allMH = NULL
allMM = NULL
allMR = NULL
allVH = NULL
allVM = NULL
allVR = NULL
#obtaining estimates for mortality, hospitalization, resident months per each quarter/year across all state
for(i in 1:numImp)
{
        currImputationH =  imputeArrVal2[,,i] #hospitalization
        currImputationM =  imputeArrVal3[,,i] #mortality
        currImputationP =  imputeArrVal1[,,i] #resident months
 #calculating rate
        diffH =  imputeArrPre1[,,i]
        diffM =  imputeArrPre2[,,i]
 #averging per quarter per year for state experiencing policy change       
        meanH = apply(diffH, 2, mean)
        varH = apply(diffH, 2, var) / length(diffH[,1])
        meanM = apply(diffM, 2, mean)
        varM = apply(diffM, 2, var) / length(diffM[,1])
        meanR = apply(currImputationP, 2, mean)
        varR = apply(currImputationP, 2, var) /  length(currImputationP[,1])
        

        allMH = cbind(allMH,meanH)
        allVH = cbind(allVH,varH)
        allMM = cbind(allMM,meanM)
        allVM = cbind(allVM,varM)        
        allMR = cbind(allMR,meanR)
        allVR = cbind(allVR,varR)    
        
}

#graph for hospitilization
pdf("Fig2HospSim.pdf")
#observed hospitalization rate
plot(1:24, forStateChange, type="l", lwd=4, ylim=c(0.04,0.1), xaxt='n',xlab = "Year (quarter)", ylab="Hospitalization Rate (per patient month)")
axis(1, at=seq(1,24,4), tick = seq(1,24,4), labels=(unique(toBedHold$year)))
abline(v=8)  policy change period

#before policy change prediction using multiple imputation formula
meanAll =  apply(allImpMB,2,mean)
varAll =  apply(allImpMB,2,var)
meanVarAll = apply(allImpVB,2,mean)
meanTotalSEAll = sqrt(meanVarAll + (1 + 1/numImp)*varAll)
upper =  meanAll + meanTotalSEAll*qnorm(0.975)
lower =  meanAll - meanTotalSEAll*qnorm(0.975)
lines(1:8, meanAll, lty=2,lwd=4,col="gray")
lines(1:8, upper, lty=3,lwd=4,col="grey50")
lines(1:8, lower, lty=3,lwd=4,col="grey50")

#after policy change prediction using MI formula 
meanAllA =  apply(allMH,1,mean)
varAllA =  apply(allMH,1,var)
meanVarAllA = apply(allVH,1,mean)
meanTotalSEAllA = sqrt(meanVarAllA + (1 + 1/numImp)*varAllA)
upperA =  meanAllA + meanTotalSEAllA*qnorm(0.975)
lowerA =  meanAllA - meanTotalSEAllA*qnorm(0.975)
polygon(c(9:24,rev(9:24)),c(upperA,rev(lowerA)),col="gray",border=NA)
lines(9:24, meanAllA, lty=2,lwd=3)
lines(1:24, forStateChange, lwd=4)#observed rate
 dev.off()
 
 
###Calculating mortality plot################################## 
#using the imputed data obtained from the MCMC sampling to calculate mortality rate
pdf("Fig2MortSim.pdf")
#observed mortality rate
plot(1:24, forStateChangeM, type="l", lwd=4, ylim=c(0.12,0.18), xaxt='n',xlab = "Year (quarter)", ylab="Mortality Rate (per patient month)")
axis(1, at=seq(1,24,4), tick = seq(1,24,4), labels=(unique(toBedHold$year)))
abline(v=8)

#before policy change prediction using multiple imputation formula
meanAllM =  apply(allImpMBM,2,mean)
varAllM =  apply(allImpMBM,2,var)
meanVarAllM = apply(allImpVBM,2,mean)
meanTotalSEAllM = sqrt(meanVarAllM + (1 + 1/numImp)*varAllM)
upperM =  meanAllM + meanTotalSEAllM*qnorm(0.975)
lowerM =  meanAllM - meanTotalSEAllM*qnorm(0.975)
lines(1:8, meanAllM, lty=2,lwd=4,col="gray")
lines(1:8, upperM, lty=3,lwd=4,col="grey50")
lines(1:8, lowerM, lty=3,lwd=4,col="grey50")

#after policy change prediction using multiple imputation formula 
meanAllAM =  apply(allMM,1,mean)
varAllAM =  apply(allMM,1,var)
meanVarAllAM = apply(allVM,1,mean)

nNumImpute = 101
nNumImp = ceiling((nNumImpute - 10)/5)
meanTotalSEAllAM = sqrt(meanVarAllAM + (1 + 1/numImp)*varAllAM)
upperAM =  meanAllAM + meanTotalSEAllAM*qnorm(0.975)
lowerAM =  meanAllAM - meanTotalSEAllAM*qnorm(0.975)
polygon(c(9:24,rev(9:24)),c(upperAM,rev(lowerAM)),col="gray",border=NA)
lines(9:24, meanAllAM, lty=2,lwd=3)
lines(1:24, forStateChangeM, lwd=4)      #observed mortality rate
 dev.off()
 
 
#############calculate resident months############################################
#graph for residents
pdf("Fig2ResSim.pdf")
#observed resident months rate
plot(1:24, forStateChangeR, type="l", lwd=4, ylim=c(240,300), xaxt='n',xlab = "Year (quarter)", ylab="Average # of Resident months")
axis(1, at=seq(1,24,4), tick = seq(1,24,4), labels=(unique(toBedHold$year)))
abline(v=8)

#before policy change prediction using multiple imputation formula
meanAllR =  apply(allImpMBR,2,mean)
varAllR =  apply(allImpMBR,2,var)
meanVarAllR = apply(allImpVBR,2,mean)
meanTotalSEAllR = sqrt(meanVarAllR + (1 + 1/numImp)*varAllR)
upperR =  meanAllR + meanTotalSEAllR*qnorm(0.975)
lowerR =  meanAllR - meanTotalSEAllR*qnorm(0.975)
lines(1:8, meanAllR, lty=2,lwd=4,col="gray")
lines(1:8, upperR, lty=3,lwd=4,col="grey50")
lines(1:8, lowerR, lty=3,lwd=4,col="grey50")

#after policy change prediction using multiple imputation formula 
meanAllAR =  apply(allMR,1,mean)
varAllAR =  apply(allMR,1,var)
meanVarAllAR = apply(allVR,1,mean)
nNumImpute = 101
nNumImp = ceiling((nNumImpute - 10)/5)
meanTotalSEAllAR = sqrt(meanVarAllAR + (1 + 1/numImp)*varAllAR)
upperAR =  meanAllAR + meanTotalSEAllAR*qnorm(0.975)
lowerAR =  meanAllAR - meanTotalSEAllAR*qnorm(0.975)
polygon(c(9:24,rev(9:24)),c(upperAR,rev(lowerAR)),col="gray",border=NA)
lines(9:24, meanAllAR, lty=2,lwd=3)
lines(1:24, forStateChangeR, lwd=4)#observed resident months rate
 dev.off()


