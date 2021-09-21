#!/bin/bash
#####
# RUN FORQS RECOMBINATION SIMULATION - one chrom only
########


### SET DEFAULTS
runID=default
pops=1
popSize=1000
founderMap=""
recMap=/home/groups/dpetrov/poolseq_pipelines/HAFpipe-line/simulations/dmel_recRates_2L-100kb.csv
gens=5
burnin=0
nloci=0
lociList=""
mDists=""
selection=0
truncation=""
report_every_nGens=1
snps=/scratch/groups/dpetrov/REFERENCE/DGRP/dgrp2.subset_all.2L.snpTable.npute
outDir="."
debug=0


### USAGE
usage()
{
    echo "
	usage: HAFpipe-sim.run_forqs.sh  
 	-id | --runID )
    -p | --pops )
	-n | --popSize )
	-f | --founderMap )
	-r | --recMap )
	-g | --gens )
    -b | --burnin )
	-nl | --nofLoci )
	-ll | --lociList )
    -md | --markerDists )
	-sc | --selectionCoeff )
    -t | --truncation )
	-s | --snps ) ## expects a file (comma or tab separated) where the first column is the snp position, the second is the reference allele, the rest are founder alleles
	-v | --report )
	-o | --out )
    -d | --debug )
	-h | --help )
"
}

### COMMAND LINE PARAMS

if [ "$1" == "" ]; then usage; exit; fi

## Parse User Parameters
while [ "$1" != "" ]; do
    case $1 in
        -id | --runID )         shift
                                runID=$1
                                ;;
        -p | --pops )      	    shift
                                pops=$1
                                ;;
        -n | --popSize )        shift
                                popSize=$1
                                ;;
	    -f | --founderMap ) 	shift
				                founderMap=$1
				                ;;
        -r | --recMap )         shift
                                recMap=$1
                                ;;
        -g | --gens )           shift
                                gens=$1
                                ;;
        -b | --burnin )         shift
                                burnin=$1
                                ;;
        -nl | --nofLoci )       shift
                                nloci=$1
                                ;;
        -ll | --lociList )      shift
                                lloci=$1
                                ;;
        -md | --markerDists )   shift
                                mDists=$1
                                ;;
        -sc | --selectionCoeff ) shift
                                selection=$1
                                ;;
        -t | --truncation )     shift
                                truncation=$1
                                ;;
        -s | --snps )         	shift
                                snps=$1
                                ;;
        -v | --report )         shift
                                report_every_nGens=$1
                                ;;
       -o | --out )        	    shift
                                outDir=$1
                                ;;
       -d | --debug )           
                                debug=1
                                ;;
       -h | --help )           usage
                                exit
                                ;;
        * )                     echo unknown flag $1 ; usage
                                exit 1
    esac
    shift
done


### MAIN

if [ ! -e $outDir/forqs.${runID} ]; then mkdir -p $outDir/forqs.${runID}; fi
cd $outDir/forqs.${runID}

if [ ! $(head -1 $snps | tr ',' '\t' | cut -f2) == "Ref" ]; then echo "error! second column of SNP table must be the ref allele"; exit 1; fi

### DEFINE FILES THAT WILL BE CREATED
###########################
seed=1
msfile=ms.txt
founderfile=founder_list.txt
configFile=config.txt


####### SET UP LOCI  LISTS 
###########################
### if only a number of selected loci is supplied, turn it into a list
if [ -e $lloci ] && [ $nloci -gt 0 ]; then lloci=$(tail -n +2 $snps | cut -f1 -d',' | shuf | head -$nloci | sort -k1,1g | tr '\n' ','); fi
### if a list of loci is supplied, count them
if [ ! -e $lloci ]; then nloci=$(echo $lloci | tr ',' ' ' | wc -w); lloci=$lloci","; fi

### if a list of marker distances is supplied, find snps appx each distance away from each locus in the list
### combine markers with selected loci for an 'all loci' list

if [ ! -z $mDists ]; then 
    mloci=$(
    echo $mDists | tr ',' '\n' | while read d; do
        echo $lloci | tr ',' '\n' | grep [0-9] | while read locus; do
            marker=$(tail -n +2 $snps | cut -f1 -d',' | awk -v locus="$locus" -v markerdist="$d" '
            function abs(x){return (x < 0) ? -x : x;}; 
            BEGIN { bestdist=1000000000 };
            { dist=abs(locus-$1); 
              if( abs(dist-markerdist) < abs(bestdist-markerdist) && ($1 != locus) ){
                bestdist=dist;marker=$1
                }
            }; 
            END{print marker}')
            echo $marker
        done
    done | sort -k1,1g |  tr '\n' ','
    )
    aloci=${lloci}${mloci}; 
else mloci=""; aloci=$lloci; fi
nmarkers=$(echo $mloci | tr ',' ' ' | wc -w); 

if [ $debug -gt 0 ]; then 
echo "
selected loci: $lloci
marker loci: $mloci
all loci: $aloci
"; 
fi

######### GET INFO FROM SNP FILE AND REC MAP
#################################
nFounders=$(head -1 $snps | tr ',' '\t' | cut -f3- | wc -w)
chromLength=$(tail -n 1 $recMap | tr ' ' '\t' | cut -f1)

###############
##  PRINT SUMMARY TO SCREEN
##############
echo "
creating forqs config file for:
 -- $pops populations of size $popSize, 
 -- each founded from a mix of $nFounders homozygous lines with genotypes specified in 
        $snps,"
if [ ! -z $founderMap ]; then 
    echo "-- using initial founder frequencies in 
        $founderMap,"  
fi

echo "
 --  recombining according to rec rates in 
        $recMap,
 -- during a $burnin generation burn-in period then $gens generations of selection,
 -- with $nloci loci under selection with strength $selection:
        $lloci
 -- and $nmarkers neutral markers:
        $mloci
 -- writing results to :
        $outDir/forqs.${runID}"


#########################
### MAKE A LIST OF FOUNDERS BASED ON A SPECIFIED SET OF FREQUENCIES OR EVEN FREQUENCIES
#########################
module load R/4.0.2
if [ ! -z $founderMap ]; then
    Rscript $HOME/scripts/simulations/write_founder_list.R $founderfile $popSize $pops $founderMap
else 
    Rscript $HOME/scripts/simulations/write_founder_list.R $founderfile $popSize $pops $nFounders
fi


#######################
## CREATE MS FILE TO SPECIFY FOUNDER GENOTYPES AT ALL TRACKED SITES
###################
if [ $nloci -gt 0 ]; then 

## create ms file to:
## assign loci to positions in genome 
## assign alleles to haplotypes
echo "" > alleles.tmp
echo $aloci | tr ',' '\n' | while read pos; do
        searchstring="^"${pos}","
        tail -n +2 $snps | grep -E $searchstring  >> alleles.tmp    
done


echo "
//
segsites: $(echo $aloci | tr ',' ' ' | wc -w)
positions: "$(cat alleles.tmp | cut -f1 |  grep [0-9] | awk -v chromLength="$chromLength" '{printf "%f ",$ii/chromLength}') > $msfile

awk -F ',' -v popSize="$popSize" -v nPops="$pops" -v nFounders="$nFounders" '
BEGIN {
    nchroms=2*popSize*nPops
}
(NR==FNR){       # collect alleles for each founder at selected sites
    pos=$1;ref=$2
    for(ii=3;ii<=NF;ii++){ if($ii==ref){allele=0}else{allele=1};hapAlleles[ii-2]=hapAlleles[ii-2]allele }
}
(NR>FNR && FNR>1){
    print hapAlleles[$3]; print hapAlleles[$3]
}' alleles.tmp $founderfile >> $msfile

fi

###############
##  SET UP CONFIG FILE 
##############

echo "
# forqs config file

# set up population
PopulationConfigGenerator_ConstantSize pcg
    population_size = $popSize
    generation_count = $(( $gens + $burnin ))
    population_count = $pops
    chromosome_pair_count = 1
    chromosome_lengths = $chromLength 
    fitness_function = qt " > $configFile 


echo "
###############
# set recomb map
###############
# generate recombination positions based on a genetic map
RecombinationPositionGenerator_RecombinationMap rpg_map
    filename = $recMap
  
  # report population haplotypes after each set of n generations  
  #  Reporter_Population reporter_population
  #  update_step = $report_every_nGens
" >> $configFile

##############
## ADD SELECTED LOCI
###############


if [ $nloci -gt 0 ]; then 


## add selected loci individually to forqs config file
cat alleles.tmp | cut -f1 -d',' |  grep [0-9] | grep -n "" | grep -f <(echo $lloci | tr ',' '\n' | grep [0-9]) | \
awk '{split($1,parts,":"); print "Locus l"parts[1]"_selected\n\tchromosome = 1\n\tposition = "parts[2]"\n"}' >> $configFile

cat alleles.tmp | cut -f1 -d',' |  grep [0-9] | grep -n "" | grep -f <(echo $mloci | tr ',' '\n' | grep [0-9]) | \
awk '{split($1,parts,":"); print "Locus l"parts[1]"_marker\n\tchromosome = 1\n\tposition = "parts[2]"\n"}' >> $configFile

echo " 
LocusList all_loci " >> $configFile
cat $configFile | grep "^Locus l" | sed 's/Locus //' | awk 'BEGIN{locistr="\tloci ="};{locistr=locistr" "$1};END{print locistr"\n"}' >> $configFile

echo " 
LocusList selected_loci " >> $configFile
cat $configFile | grep "^Locus l.*_selected$" | sed 's/Locus //' | awk 'BEGIN{locistr="\tloci ="};{locistr=locistr" "$1};END{print locistr"\n"}' >> $configFile

echo " 
LocusList marker_loci " >> $configFile
cat $configFile | grep "^Locus l.*_marker$" | sed 's/Locus //' | awk 'BEGIN{locistr="\tloci ="};{locistr=locistr" "$1};END{print locistr"\n"}' >> $configFile

echo " 
VariantIndicator_File variants
    msfile = $msfile
    loci = all_loci
" >> $configFile


## set selection coefficients for each selected locus individually
if [ ! -z "$burnin" ]; then
    echo "QuantitativeTrait_IndependentLoci qt_selection" >> $configFile
else
    echo "QuantitativeTrait_IndependentLoci qt" >> $configFile
fi

cat $configFile | grep "^Locus l.*_selected$" | sed 's/Locus //' | awk -v sc="$selection" '
{print "\tqtl = "$1" 1 "1+sc/2" "1+sc}
' >> $configFile
echo "    environmental_variance = 0.05
" >> $configFile

## add burnin period and generate composite qt trajectory
if [ ! -z "$burnin" ]; then
    echo "QuantitativeTrait_IndependentLoci qt_burnin" >> $configFile
    cat $configFile | grep "^Locus l.*_selected$" | sed 's/Locus //' | awk '{print "\tqtl = "$1" 1 1 1"}' >> $configFile
    echo "
    " >> $configFile

    echo "QuantitativeTrait_GenerationComposite qt
    generation:quantitative_trait = 0 qt_burnin
    generation:quantitative_trait = 3 qt_selection
    " >> $configFile

fi

## create fitness function
if [ ! -z "$truncation" ]; then
echo "FitnessFunction_TruncationSelection fitness
    quantitative_trait = qt
    proportion_selected = $truncation
    " >> $configFile
fi

#report allele frequencies at all tracked sites
echo "Reporter_AlleleFrequencies reporter_allele_freqs
    locus_list = all_loci

"  >> $configFile

if [ $nloci -eq 1 ]; then
    selected_ix=$(cat alleles.tmp | grep ',' | grep -n $lloci | cut -f1 -d":")
	af0=$(tail -n +5 $msfile | awk -F '' -v selix="$selected_ix" '{sum=sum+$selix}END{print sum/NR}')
	w1=$( awk -v sc="$selection" 'BEGIN{print sc/2 + 1}' )
	w2=$( awk -v sc="$selection" 'BEGIN{print sc + 1}' )
	echo "Reporter_DeterministicTrajectories reporter_deterministic_trajectories
    initial_allele_frequency = $af0
    w0 = 1
    w1 = $w1
    w2 = $w2
	" >> $configFile
fi

fi





############
## PUT IT ALL TOGETHER
###########
# SimulatorConfig is the top-level module, and must appear last in the
# configuration file.  This module specifies which of the previously
# specified primary modules are plugged into the main simulator.

echo "
SimulatorConfig
    output_directory = $outDir/forqs.${runID}/results
    population_config_generator = pcg
    recombination_position_generator = rpg_map
#    reporter = reporter_population
    seed = $seed"  >> $configFile
if [ $nloci -gt 0 ]; then echo "    quantitative_trait = qt" >> $configFile; fi
if [ ! -z $truncation ]; then echo "    quantitative_trait = fitness" >> $configFile; fi
if [ $nloci -gt 0 ]; then echo "    reporter = reporter_allele_freqs
    variant_indicator = variants" >> $configFile; 
	if [ $nloci -eq 1 ]; then 
	echo "	reporter = reporter_deterministic_trajectories" >> $configFile	
	fi
fi

if [ $debug -gt 0 ]; then 
    echo "setup written to $configFile - run forqs with:"
    echo "forqs $outDir/forqs.${runID}/$configFile > $outDir/forqs.${runID}/log"
else
    echo "setup written to $configFile - now running forqs..."
    if [ -e $outDir/forqs.${runID}/results ]; then rm -r $outDir/forqs.${runID}/results; fi
    forqs $configFile > log

    echo "--> results written to 
    $outDir/forqs.${runID}
    "

    echo $lloci | tr ',' '\n' | grep [0-9] | while read locus; do 
        mv results/allele_frequencies_chr1_pos${locus}.txt \
        results/af_selected_chr1_pos${locus}.txt; 
    done
    echo $mloci | tr ',' '\n' | grep [0-9] | while read locus; do 
        mv results/allele_frequencies_chr1_pos${locus}.txt \
        results/af_marker_chr1_pos${locus}.txt; 
    done
fi
