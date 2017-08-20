#######################Synthetic Control with Simulated Data##################
##############################################################################
require(doBy)
require(Synth)

library(doBy)
library(Synth)
#loading simulated data
dataUse = read.csv("simulatedDataSetForPaper.csv",header=TRUE,row.names=1)
numBefore = 8
numAfter = 16

# Split states by bedhold policy
toBedHold = subset(dataUse,dataUse$BedHoldChange==1)
noBedHold = subset(dataUse,dataUse$BedHoldChange==0)

#switching to rate log scale
toBedHoldF = toBedHold
toBedHoldF[,1] = (toBedHoldF[,1]/toBedHoldF[,3])
toBedHoldF[,2] = (toBedHoldF[,2]/toBedHoldF[,3])
toBedHoldF[,3] = (toBedHoldF[,3])

#creating quarter sequence to use with the synthetic method function
tofaclength2 = summaryBy(quarter ~ facid, data = toBedHold, FUN = length)[,2]
tofaclen = NULL
for(i in 1:length(tofaclength2)) tofaclen = c(tofaclen, seq(1:tofaclength2[i]))
toBedHoldF = cbind(toBedHoldF, tofaclen)
colnames(toBedHoldF)[length(toBedHoldF[1,])] = "nofaclen"
# noBedHoldF
#switching to rate log scale for no bed hold facilities
noBedHoldF = noBedHold
noBedHoldF[,1] = (noBedHold[,1]/noBedHold[,3]) #acute hospitalization
noBedHoldF[,2] = (noBedHold[,2]/noBedHold[,3]) #death
noBedHoldF[,3] = (noBedHold[,3])                         #person months


nofaclength2 = summaryBy(quarter ~ facid, data = noBedHoldF, FUN = length)[,2]
nofaclen = NULL
for(i in 1:length(nofaclength2)) nofaclen = c(nofaclen, seq(1:nofaclength2[i]))
noBedHoldF = cbind(noBedHoldF, nofaclen)

###Start preparing for synthetic
#Used to running the script in parallel there are 413 facilities in the state that changed policy
#running them 50 at a time using minFacID the lower number and maxFacID the upper number
args=commandArgs()
if(length(args)==0)
{    
	print("No arguments supplied.");    
	minFacID = 1;    
	maxFacID = 50;
}else
{    
	for(lenJ in 9:length(args))
	{         
		eval(parse(text=args[[lenJ]]))
		    
	}
}
#file to export the results
fileName = sprintf("simSynthD_%d_%d.csv",minFacID,maxFacID) 
#total number of controls to use for synthetic control for each oneof the 413 facilities
numControls = 40

#removing certain facilities that make synthetic control program crash
if(minFacID == 101)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid %in% c(951,1547)),]
}
#if(minFacID == 151)
#{
#    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid == 1097),]
#}
if(minFacID == 251)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid %in%  c(925)),]
}
if(minFacID == 301)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid %in% c(2308,1862,763)),]
}
if(minFacID == 51)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid %in%  c(1243)),]
}
if(minFacID == 1)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid %in% c(1865,951)),]
}

if(minFacID == 401)
{
    noBedHoldF = noBedHoldF[-which(noBedHoldF$facid == 615),]
}

#end removal of units
#identifying the periods before any policy change occcured in the states that didn't change policy
lessThan8no = which(noBedHoldF$nofaclen <=8)
sizeNoBedhold = tapply(noBedHoldF$nhmod[lessThan8no], noBedHoldF$facid[lessThan8no],mean)

#identifying the periods before any policy change occcured in the states that changed policy
lessThan8to = which(toBedHoldF$nofaclen <=8)
sizeToBedhold = tapply(toBedHoldF$nhmod[lessThan8to], toBedHoldF$facid[lessThan8to],mean)


facidUse = unique(toBedHoldF$facid)
allGapMat = NULL
#for each facility that changed policy
for(j in minFacID:min(maxFacID,length(facidUse)))
{
#prepare dataset to be used by synthetic control function
  sizeVal =  sizeToBedhold[as.numeric(names(sizeToBedhold))== facidUse[j]]
  useNoBed = as.numeric(names(sizeNoBedhold[order(abs(sizeNoBedhold-sizeVal))])[1:numControls])
  newFile = rbind(toBedHoldF[toBedHoldF$facid == facidUse[j],],noBedHoldF[which(noBedHoldF$facid %in% useNoBed),])
  newFile = cbind(newFile,newFile$facid)
  colnames(newFile)[12] = "facidN"
  newFile = data.frame(newFile)
  newFile$facidN = as.character(newFile$facidN)
   print(sprintf("*********** %d ***********",facidUse[j]))
#  if(facidUse[j] ==3175)
 #{
	#next;
#}
 #running synthetic control for facility j
  dataprep.out <- dataprep(
               foo = newFile,
               predictors = c("nhmod","acute_qtr"),
               predictors.op = "mean",
               time.predictors.prior = 1:numBefore,
               special.predictors = list( list("died_qtr", 1:numBefore , "mean")),
               dependent = "died_qtr",
               unit.variable = "facid",
               unit.names.variable = "facidN",
               time.variable = "nofaclen",
               treatment.identifier = facidUse[j],
               controls.identifier = useNoBed,
               time.optimize.ssr = 1:numBefore,
               time.plot = 1:(numBefore+numAfter))

  synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")
  gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w)#calculating different from sythetic
  allGapMat = rbind(allGapMat, c(facidUse[j],1,gaps)) #adding treated facility synthetic

  #this would be used to create the "placebo test" it calculates the difference for each of the controls from the other chosen controls
  for(i in 1:numControls)
  {
    print(sprintf("###################### %d #############",useNoBed[i]))
   dataprep.out <- dataprep(
   foo = newFile[newFile$facid != facidUse[j],],
   predictors = c("nhmod","acute_qtr"),
   predictors.op = "mean",
   time.predictors.prior = 1:numBefore,
   special.predictors = list( list("died_qtr", 1:numBefore , "mean")),
   dependent = "died_qtr",
   unit.variable = "facid",
   unit.names.variable = "facidN",
   time.variable = "nofaclen",
   treatment.identifier = useNoBed[i],
   controls.identifier = useNoBed[-i],
   time.optimize.ssr = 1:numBefore,
   time.plot = 1:(numBefore+numAfter))

   synth.out <- synth(data.prep.obj = dataprep.out, method = "BFGS")
   gaps <- dataprep.out$Y1plot - (dataprep.out$Y0plot %*% synth.out$solution.w) #calculating difference from synthetic
   allGapMat = rbind(allGapMat, c(facidUse[j],0,gaps)) #adding control facility synthetic
  }
  #wrtitng results to a file
  write.csv(allGapMat,file=fileName)
}
