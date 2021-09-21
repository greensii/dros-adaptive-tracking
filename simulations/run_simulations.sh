###########
## simulations for orchard 2014 paper
#############

forqsScript=forqs_wrapper.sh


#### goals
#############

# simulate 10 replication populations
pops=10

# each with a population of 100k
n=100000

# from the founders used in the cages
snps=inbredv2_withHets.subset_orch14.2L.snpTable.npute

# start with even founder frequencies
founderMap=""

# add a 3-generation burn-in period (ie. neutral recombination; no selection)
burnin=3

# set selection coefficients from .01 to .5
selCoefs=( 0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 ) 

# simulate selection for 5 generations
gens=5

# simulate truncation selection with top 25% of the pop reproducing
trunc=.25

######################
## run forqs

if [! -d results ]; then mkdir results; done


site_groups=( 100 50 20 10 ) 

for ngroups in ${site_groups[*]}; do

	site_list=sim_sites.${ngroups}_groups.csv
	cat sim_sites.csv | shuf | awk -v ngroups="$ngroups" '{g=NR % ngroups; groups[g]=groups[g]""$1"\t" }; END{for(g in groups){print groups[g]}}' > $site_list

	cat $site_list | while read line; do
		af0=$(echo $line | awk '{for(ii=1;ii<=NF;ii++){split($ii,parts,":");af0=af0""parts[1]"."};print substr(af0,1,length(af0)-1)}')
		pos=$(echo $line | awk '{for(ii=1;ii<=NF;ii++){split($ii,parts,":");split(parts[2],parts,",");pos=pos""parts[1]","};print substr(pos,1,length(pos)-1)}')
	

		for sel in ${selCoefs[*]}; do
	
			cmd="$forqsScript \
	-id 2L_af0-${af0}_sel-${sel}_trunc-${trunc} \
	-p $pops \
	-n $n \
	-b $burnin \
	-g $gens \
	-ll $pos \
	-sc $sel \
	-md 5000 \
	-s $snps \
	-t $trunc \
	-o results "

			echo $af0","$sel

			#	echo $cmd
			eval $cmd
	
		done
	done
done
