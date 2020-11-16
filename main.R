# MAIN

## locations
#localDir="~/Documents/cage_data/orchard_2014";
localDir="/scratch/users/greensi/orchard_2014"
setwd(localDir)

teamdriveDir="teamdrive:Dmel_cageExperiments_SGreenblum/2014-orchard-dmel"
#scriptsDir="~/Documents/GitHub//dros-adaptive-tracking"
scriptsDir="~/scripts/GitHub/dros-adaptive-tracking"
looDir="Rdata//glm_loo"

## files
snpFile="snps/inbredv2_withHets.subset_orch14.CHROM.snpTable.numeric"
fstFile="fst_perChrom.orch14_Ecages_baseline.Rdata"
HAFsFile="Rdata//HAFs.orch14_Ecages.Rdata"
baselineFile="Rdata/ds5xHAFs.orch14_baseline.Rdata"
glmFile="Rdata//glm.Ecages.Rdata"


## globals
chroms=c("2L","2R","3L","3R","X")
comparisons=c("2_1","3_2","4_3","5_4","5_1")
timesegs=c("Timepoint 1 --> 2", "Timepoint 2 --> 3", "Timepoint 3 --> 4","Timepoint 4 --> 5","Timepoint 1 --> 5")
colorScheme=c("BH-FDR<.2"="lightgray",
              "BH-FDR<.05,\neffect-size>2%"="coral",
              "empirical-FDR<.05,\neffect-size>2%"="red")

## parameters
fdrThreshs=c(.2,.05,.01) ## maximum fdr-corrected pvalue for difference in allele frequency between treatments for a site to be considered significantly diverged
esThreshs=c(0,0.02,0.02) ## minimum mean allele frequency difference between treatments for a site to be considered significantly diverged (combined with p-value)
windowSize=500
windowShift=100
maxClusterGap=100
linkageThresh=0.03
maxSNPpairDist=3000000
poolSize=100
maxthreads=20 ## max number of parallel processes to run (actual num will not exceed available cores detected below)

## FUNCTIONS
setwd(scriptsDir)
#system("git pull",intern=T)
source("workflow_functions.R")
source("plotting_functions.R")
source("helper_functions.R")
source("load_packages.R")
setwd(localDir)

results=RunFullWorkflow(HAFsFile,glmFile,snpFile,comparisons,1:10,fdrThreshs,esThreshs,windowSize,windowShift,maxClusterGap,maxSNPpairDist,linkageThresh,maxthreads)