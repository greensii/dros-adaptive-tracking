# dros-adaptive-tracking

### genomic analyses for experimental evolution of 10 replicate drosophila melanogaster populations adapting to seasonal shifts

_Seth M. Rudman, Subhash Rajpurohit, Sharon Greenblum, Ozan Kiratli, Martin M. Turcotte, Dmitri A. Petrov, Jonathan Levine, Paul Schmidt_

#### Experimental Design:
+ In June of 2014, large replicate cages located in a Pennsylvania orchard were each seeded with 500 male and 500 female D melanogaster from a single reconstructed outbred founder population. 
  + The founder population was constructed from ~80 homozygous inbred lines collected in spring in PA.
  + 4 biological replicate samples of 100 individuals were collected from the founder population at the time of cage seeding and sequenced at ~50x coverage.
+ 10 replicate cages were allowed to freely evolve.
  + Each replicate population was allowed to evolve for ~15 generations
  + At 5 timepoints between June and November 100 female individuals from each cage were sampled for pooled sequencing 
+ Whole genome pool-seq was performed at ~5x coverage per sample, and then haplotype inference was used to calculate allele frequencies (Tilk et al, 2019).

#### Summary of results:
+ Rapid adaptation is strong enough to prompt detectable genomic shifts within a population on monthly timescales
  + average Fst vs baseline shows divergence is greater than expected variation due to sampling or technical error, and changes over time
+ Replicate populations in the same environment shift in parallel
  + leave-one-out GLM analysis shows that allele frequency shifts in 9 cages are predictive of 10th left-out cage 
  + full GLM analysis shows that per-site parallelism is significantly stronger than empirical null (from shuffling sample timepoint labels)
+ The parallel response is polygenic (or at least oligogenic)
  + Significantly parallel sites are found on every chromosome
  + Significant sites are over-dispersed, consistent with selection+linkage
  + After accounting for linkage, close to 150 independent clusters of linked sites shift significantly 
  + Most clusters are not linked to inversions
+ Selection fluctuates on sub-seasonal timescales
  + LOO GLM analysis shows direction of allele frequency shift in 10th left-out cage often flips between time segments
  + Few regions are enriched for parallelism over multiple time segments
  + SNPs in the regions that do span multiple time segments often switch direction
  
#### Workflow
workflow steps performed by the functions in this repo:  
+  _get_glm_pvals_: fit allele freqs at each SNP to GLM and extract p-values for significance of sampling timepoint  
+  _get_af_shifts_: calculate average allele frequency shift for each cage over each time interval
+  _get_sig_sites_: identify sites that meet FDR and effect size thresholds  
+  _score_wins_: score sliding genomic windows for enrichment for significantlhy parallel sites  
+  _get_win_fdr_: compare distribution of window scores to empirical null, identify significantly enriched windows  
+  _cluster_wins_: merge overlapping enriched windows  
+  _associate_snps_to_clusters_: identify the significant SNPs in each cluster  
+  _find_snp_pairs_: create pairs of significant SNPs, find matched control pairs  
+  _calc_Rsq_for_snp_pairs_: calculate linkage (squared correlation coef) for selected SNP pairs  
+  _merge_linked_clusters_: merge consecutive clusters with high SNP-pair linkage  
  
Note that this same workflow was performed using all 10 cages, and also in a 10-fold leave-one-cage-out cross-validation

#### Input Files:
Example input files for the first 1,000 sites on chromosome 2L are included in the test_data directory:
+  _founder_genotypes_1ksites.2L.snpTable.numeric_ is a comma-separated text file specifying the genotype of each founder line (columns) at each segregating site (rows). The first column lists the position of the site, the rest of the columns code for genotypes as follows:
	+   0 is homozygous reference
	+   1 is homozygous alternate
	+   0.5 is heterozygous (should be very few of these with inbred lines)
	+   -1 is missing genotype call 
+  _HAFs_1ksites.Rdata_: contains the variables 'sites', 'samps', and 'afmat'
	+   sites is a dataframe corresponding to the rows of afmat, with columns 'chrom' and 'pos'
	+   samps is a dataframe corresponding to the columns of afmat, with columns 'sampID', 'tpt', and 'cage'
	+   afmat is a numeric matrix of haplotype-inferred allele frequencies for each site (rows) in each sample (columns)
+  _glm_1ksites.Rdata_: contains the variable 'df.glm', a dataframe with GLM coefficients and p-values corresponding to the significance of parallelism at each time segment for each tested site

Author: Sharon Greenblum, 2018-2019 Stanford University  
To cite, please reference [![DOI](https://zenodo.org/badge/311184324.svg)](https://zenodo.org/badge/latestdoi/311184324)
