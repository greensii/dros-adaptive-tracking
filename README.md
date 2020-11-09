# dros-adaptive-tracking

### genomic analyses for experimental evolution of 10 replicate drosophila populations adapting to seasonal shifts

_Seth M. Rudman, Subhash Rajpurohit, Sharon Greenblum, Ozan Kiratli, Martin M. Turcotte, Dmitri A. Petrov, Jonathan Levine, Paul Schmidt_

#### Experimental Design:
+ In June of 2014, large replicate cages located in a Pennsylvania orchard were each seeded with 500 male and 500 female D melanogaster from a single reconstructed outbred founder population. 
  + The founder population was constructed from ~80 homozygous inbred lines collected in spring in PA.
  + 4 biological replicate samples of 100 individuals were collected from the founder population at the time of cage seeding and sequenced at ~50x coverage.
+ 10 replicate cages were allowed to freely evolve.
  + Each replicate population was allowed to evolve for ~15 generations
  + At 5 timepoints between June and November 100 female individuals from each cage were sampled for pooled sequencing 
+ Whole genome pool-seq was performed at ~5x coverage per sample, and then haplotype inference was used to calculate allele frequencies.

#### Summary of results:
+ Rapid adaptation is strong enough to prompt detectable genomic shifts within a population on monthly timescales
  + PCA and plots and average Fst vs baseline show divergence increasing over time
+ Replicate populations in the same environment shift in parallel
  + leave-one-out GLM analysis shows that allele frequency shifts in 9 cages are predictive of 10th left-out cage 
  + GLM analysis shows that per-site parallelism is significantly stronger than empirical null (from shuffling sample timepoint labels)
+ The parallel response is polygenic (or at least oligogenic)
  + Significantly parallel sites are found on every chromosome
  + Thousands of 500-SNP (~40kb) genomic windows are enriched for parallelism
  + Enriched regions are both inside and outside inversions, and not predominantly near centromeres where recomb is suppressed
  + After accounting for linkage, close to 100 unique independent site groups shift significantly 
+ Selection fluctuates on sub-seasonal timescales
  + LOO GLM analysis shows direction of allele frequency shift in 10th left-out cage often flips between time segments
  + Few regions are enriched for parallelism over multiple time segments
  + SNPs in the regions that do span multiple time segments often switch direction
  + (PCA plots suggest genome-wide reversion towards baseline (spring-adapted) state at the end of the sampling season?). 
  
Author: Sharon Greenblum, 2018-2019 Stanford University
