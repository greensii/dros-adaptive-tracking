## files
snpFile="data/snps/inbredv2_withHets.subset_orch14.CHROM.snpTable.numeric"
fstFile="data/fst_perChrom.orch14_Ecages_baseline.Rdata"
HAFsFile="data/HAFs.orch14_Ecages.Rdata"
glmFile="data/glm.Ecages.Rdata"
baselineFile="data/ds5xHAFs.orch14_baseline.Rdata"
invFile="data/inversions/dmel_inversion_breakpoints.fromCorbittDetigandHartl2012.txt"

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

