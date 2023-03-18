# MAIN

## setup - edit this to your desired paths
localDir= ### ie. "~/data/orchard_2014";
scriptsDir= ### ie. "~/GitHub/dros-adaptive-tracking"

## FUNCTIONS
setwd(scriptsDir)
source("workflow_functions.R")
source("plotting_functions.R")
source("helper_functions.R")
source("load_packages.R")
setwd(localDir)

## PARAMETERS
source(configFile)

### RUN WORKFLOW

## for each leave-one-cage-out round
results.loo=lapply(cages,function(dropCage){
  loofile=paste0("glm_loo",cages[dropCage],".orch14.Rdata")
  RunFullWorkflow(HAFsFile,loofile,snpFile,comparisons,cages[-dropCage],fdrThreshs,esThreshs,windowSize,windowShift,maxClusterGap,maxSNPpairDist,linkageThresh,maxthreads)
})
list[shifts.loo,medians.loo]<-get_loo_shifts(results.loo,HAFsFile,snpFile,comparisons,timesegs)
save(results.loo,shifts.loo,medians.loo,file="Rdata/results_loo.orch14.Rdata")

## for all cages together
results=RunFullWorkflow(HAFsFile,glmFile,snpFile,comparisons,cages,fdrThreshs,esThreshs,windowSize,windowShift,maxClusterGap,maxSNPpairDist,linkageThresh,maxthreads)
save(results,file="Rdata/results.orch14.Rdata")
