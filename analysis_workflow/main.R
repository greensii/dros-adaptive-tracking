# MAIN

## setup - local
localDir="~/Documents/cage_data/orchard_2014";
scriptsDir="~/Documents/GitHub//dros-adaptive-tracking"

## setup - cluster
localDir="/scratch/users/greensi/orchard_2014"
scriptsDir="~/scripts/GitHub/dros-adaptive-tracking"

## FUNCTIONS
setwd(scriptsDir)
#system("git pull",intern=T)
source("workflow_functions.R")
source("plotting_functions.R")
source("helper_functions.R")
source("load_packages.R")
setwd(localDir)

## PARAMETERS
source(configFile)

## globals
chroms=c("2L","2R","3L","3R","X")
cages=1:10
comparisons=c("2_1","3_2","4_3","5_4","5_1")
timesegs=c("Timepoint 1 --> 2", "Timepoint 2 --> 3", "Timepoint 3 --> 4","Timepoint 4 --> 5","Timepoint 1 --> 5")
colorScheme=c("lightgray","coral","red")
names(colorScheme)= c(paste0("BH-FDR<",fdrThreshs[1]),
                      paste0("BH-FDR<",fdrThreshs[2],",\neffect-size>",esThreshs[2]*100,"%"),
                      paste0("BH-FDR<",fdrThreshs[3],",\neffect-size>",esThreshs[3]*100,"%")
)

### RUN WORKFLOW

## for each leave-one-cage-out round
results.loo=lapply(cages,function(dropCage){
  loofile=paste0("data/glm_loo/glm.Ecages",paste0(sort(as.character(cages)[-dropCage]),collapse=""),".Rdata")
  RunFullWorkflow(HAFsFile,glmFile,snpFile,comparisons,cages,fdrThreshs,esThreshs,windowSize,windowShift,maxClusterGap,maxSNPpairDist,linkageThresh,maxthreads)
})
list[shifts.loo,medians.loo]<-get_loo_shifts(results.loo,HAFsFile,snpFile,comparisons,timesegs)
save(results.loo,shifts.loo,medians.loo,file="Rdata/results_loo.orch14.Rdata")

## for all cages together
results=RunFullWorkflow(HAFsFile,glmFile,snpFile,comparisons,cages,fdrThreshs,esThreshs,windowSize,windowShift,maxClusterGap,maxSNPpairDist,linkageThresh,maxthreads)
save(results,file="Rdata/results.orch14.Rdata")



## backup to teamdrive
backupDir="schmidt-petrov-td:data/orchard_2014"
system("rclone copy . ",backupDir,"/notebooks --include *.ipynb --max-depth 1",intern=T)
system(paste0("rclone copy data ",backupDir,"/Rdata --exclude OLD/"),intern=T)
