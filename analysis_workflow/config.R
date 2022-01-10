## files
snpFile="test_data/founder_lines_1ksites.2L.snpTable.numeric"
HAFsFile="test_data/HAFs_1ksites.Rdata"
glmFile="test_data/glm_1ksites.Rdata"

## globals
chroms=c("2L","2R","3L","3R","X")
cages=1:10
comparisons=c("2_1","3_2","4_3","5_4","5_1")
timesegs=c("Timepoint 1 --> 2", "Timepoint 2 --> 3", "Timepoint 3 --> 4","Timepoint 4 --> 5","Timepoint 1 --> 5")
colorScheme=c("BH-FDR<.2"="lightgray",
              "BH-FDR<.05,\neffect-size>2%"="coral",
              "BH-FDR<.01,\neffect-size>2%"="red")

## parameters
fdrThreshs=c(.2,.05,.01) ## maximum fdr-corrected pvalue for difference in allele frequency between treatments for a site to be considered significantly diverged
esThreshs=c(0,0.02,0.02) ## minimum mean allele frequency difference between treatments for a site to be considered significantly diverged (combined with p-value)
windowSize=500
windowShift=100
maxClusterGap=100
linkageThresh=0.03
maxSNPpairDist=3000000
poolSize=100

