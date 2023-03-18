source("load_packages.R")

### TODO
## add way to run on subset of samples (ie only samples from a specific treatment group)? (this is separate from dropping a single sample/cage)
## add option to shuffle mainEffect labels within reps

## FUNCTIONS
#################
parse_cl_args=function(args){
  # Create a parser
  p <- arg_parser("Calculate GLM coefficients and p-values")
  
  # Add command line arguments
  p <- add_argument(p, "HAFs", help=".Rdata file containing 3 objects: afmat(matrix), samples(data.frame), sites(data.frame)")
  p <- add_argument(p, "--readDepth", help=".RDS file containing a matrix with same dimensions as afmat in HAFs, giving the raw read depth per site/sample", default=NA)
  p <- add_argument(p, "--effectiveCov", help="either a single number to be used as the effective coverage for every site/sample, or \na .csv file with column names ['sampID','chrom','ec'] containing an estimate of effective coverage per chrom/sample", default=NA)
  p <- add_argument(p, "--dropRep", help="ID of replicate to drop (when running leave-one-out)", default=NA)
  p <- add_argument(p, "--poolSize", help="number individuals sampled per pool", default=100,type="integer")
  p <- add_argument(p, "-- mainEffect", help="calculate p-values for all pairwise comparisons of groups in this sample metadata column", default="tpt")
  p <- add_argument(p, "--repName", help="name of the column in the sample metadata table that identifies replicate IDs", default="cage")
  p <- add_argument(p, "--testNsites", help="run GLM on a random subset of N sites", default=NA,type="integer")
  p <- add_argument(p, "--ncores", help="run GLM in parallel using mclapply with this many cores", default=1,type="integer")
  p <- add_argument(p, "--outDir", help="write all results to this directory; will be created if it doesnt already exist", default=".")
  
  # Parse the command line arguments
  args <- parse_args(p)
  
  ## check HAFs
  if(! file.exists(args$HAFs) ){cat("HAFs file",args$HAFs,"does not exist\n***EXITING***\n");quit() }
  if(! grepl("\\.Rdata$",args$HAFs) ){cat("HAFs file",args$HAFs,"must be saved as an .Rdata file\n***EXITING***\n");quit() }
  
  ## must supply either readDepth or effectiveCov; effCov can be a file or a single integer
  ## readDepth should be a file containing a matrix with the same dimensions as afmat in HAFs, giving the raw read depth per site/sample
  ## if effectiveCov is a number, it will be used as the effective coverage for every site/sample
  ## if effectiveCov is a file, it should contain an estimate of effective coverage per chrom/sample
  if(is.na(args$readDepth)){
    if(is.na(args$effectiveCov)){
      cat("either --readDepth or --effectiveCov must be supplied.\n***EXITING***\n");quit()
    } else{
      if(is.na(as.numeric(args$effectiveCov)) & !file.exists(args$effectiveCov)){
        cat("--effectiveCov must either be a single integer or a file\n***EXITING***\n");quit()
      } else {
        if(!file.exists(args$effectiveCov)){args$effectiveCov=as.integer(args$effectiveCov)}
      }
    }
  } else {
    if(!file.exists(args$readDepth)){cat("readDepth file",args$readDepth,"does not exist\n***EXITING***\n");quit()}
  }
  if(!dir.exists(args$outDir)){cat("making outDir",args$outDir,"\n");dir.create(args$outDir,recursive = TRUE)}
  
  ## check ncores,testNsites,poolSize
  if (is.na(as.numeric(args$ncores))) {cat("ncores must be an integer\n***EXITING***\n");quit()}
  if (!is.na(args$testNsites) && is.na(as.numeric(args$testNsites))) {cat("testNsites must be an integer\n***EXITING***\n");quit()}
  if (is.na(args$poolSize)) {cat("poolSize must be an integer\n***EXITING***\n");quit()}
  
  return(args)
}
set_up_sampData=function(sampIX,samps,model.vars,cmpAll){
  sampData=samps[sampIX,]%>%dplyr::select(all_of(model.vars));
  varCols=which(apply(sampData,2,function(x){length(unique(x))>1}))
  sampData=sampData%>%dplyr::select(all_of(varCols))
  if(!is.null(cmpAll)){sampData[[cmpAll]]=factor(sampData[[cmpAll]])}
  return(sampData)
} 
set_up_depthData=function(args,sites,samps){
  if(!is.na(args$readDepth)){
    cat("using raw read depth for binomial counts\n")
    rd=readRDS(args$readDepth)
    ss=match(samps$sampID,colnames(rd))
    if(sum(is.na(ss)>0)){cat("--readDepth file is missing samples:",unique(samps$sampID[is.na(ss)]),"\n***EXITING***\n");quit()}
    rd=rd[,ss]
    if(nrow(rd)!=nrow(sites)){cat("--readDepth file must have the same number of rows as sites/afmat in --HAFs file\n***EXITING***\n");quit()}
    ### NOTE: not checking whether sites in rd and sites actually match, just that nof rows are the same !
  } else{
    if(!is.numeric(args$effectiveCov)){
      cat("using per-chromosome effective coverage for binomial counts\n")
      df.ec=fread(args$effectiveCov);
      if(length(intersect(c("sampID","chrom","ec"),colnames(df.ec)))<3){
        cat("--effectiveCov file must have the following column names: sampID,chrom,ec\n***EXITING***\n");quit()
      }
      df.ec$chrom <- factor(df.ec$chrom)
      mat.ec = t(xtabs(ec ~ sampID + chrom, df.ec))
      cc=match(sites$chrom,rownames(mat.ec))
      ss=match(samps$sampID,colnames(mat.ec))
      if(sum(is.na(cc))>0){cat("--effectiveCov file is missing chromosomes:",unique(sites$chrom[is.na(cc)]),"\n***EXITING***\n");quit()}
      if(sum(is.na(ss))>0){cat("--effectiveCov file is missing samples:",unique(samps$sampID[is.na(ss)]),"\n***EXITING***\n");quit()}
      rd=round(mat.ec[cc,ss])
    } else{
      cat("using constant value for binomial counts\n")
      rd=matrix(args$effectiveCov,nrow(sites),nrow(samps))
    }
  }
  return(rd)
}
extract_coef_pval=function(model,cmpAll=NULL,dontReport=NULL){
  if(!is.null(cmpAll)){
    model.multcomp=summary(eval(parse(text=paste0("glht(model, mcp(",cmpAll,"='Tukey'))")))) 
    cp=cbind(coefficients(model.multcomp),model.multcomp$test$pvalues)
    row.names(cp)=gsub(" - ","_",row.names(cp))
  } else{cp=summary(model)$coefficients[-1,c(1,4),drop=FALSE];}
  if(!is.null(dontReport)){cp=cp[grep(paste0("(",paste0(dontReport,collapse="|"),")"),row.names(cp),invert = TRUE),,drop=FALSE];}
  return(cp)
}
fit_GLM_one = function(af.site,rd.site,sampData,formulaString,poolSize,cmpAll=NULL,dontReport=NULL){
  Neff=((poolSize*2*rd.site)-1)/(poolSize*2+rd.site);
  cts=cbind(round(Neff*af.site),round(Neff*(1-af.site)))
  model=glm(as.formula(formulaString),family="quasibinomial",data=sampData) 
  cp=extract_coef_pval(model,cmpAll,dontReport)      
}
fit_GLM_all=function(siteIX,sampIX,samps,model.vars,poolSize,cmpAll=NULL,dontReport=NULL){
  
  sampData=set_up_sampData(sampIX,samps,model.vars,cmpAll)
  
  formulaString=paste0(colnames(sampData),collapse=" + ")
  formulaString=paste0("cts ~ ",formulaString," + 1")
  
  cat("Model Formula is: \n",formulaString,"\n")
  nSites=length(siteIX);cat("there are ",nSites," sites\n")
  do.call(rbind,mclapply(siteIX,function(ix){
    
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    
    cp=fit_GLM_one(afmat[ix,sampIX],rd[ix,sampIX],sampData,formulaString,poolSize,cmpAll,dontReport)
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.",row.names(cp)),paste0("p.",row.names(cp)))
    return(results)},mc.cores=args$ncores))
  }

##########
## MAIN
args <- suppressWarnings(parse_cl_args(commandArgs(trailingOnly=TRUE)))
registerDoMC(cores=args$ncores)
cat("RUNNING calc_GLM.R with the following parameters:\n")
print(args);

# load af data
cat("loading HAFs\n")
load(args$HAFs)  ## should contain sites, samps, afmat

## get chroms and sites
chroms=unique(sites$chrom)
if(!is.na(args$testNsites)){
  siteIX=sample(1:nrow(sites),args$testNsites)
} else {siteIX=1:nrow(sites)}

# set up samp info
samps$cage=samps[[args$repName]]
samps$cage=factor(samps$cage)
cage_set=sort(levels(samps$cage)) 
if(!is.na(args$dropRep)){cage_set=cage_set[cage_set!=args$dropRep]}
sampIX=which(samps$cage %in% cage_set);

# set up read depth // effective coverage info
rd=set_up_depthData(args,sites,samps)

# calc GLM
cat("calculating GLM using main effect:",args$mainEffect,
    "\nusing",length(sampIX),"samples from replicates:",cage_set,"\n")
system.time({
    df.glm=cbind(sites[siteIX,],fit_GLM_all(siteIX,sampIX,samps,model.vars=c(args$mainEffect,args$repName),
                                            poolSize=args$poolSize,cmpAll=args$mainEffect,dontReport=args$repName))
})
cat("finished\n")

df.glm <- df.glm %>% mutate_at(vars(starts_with("p.")),function(x){x[x==0]=NA;return(x)}) 
filename=paste0(args$outDir,"/glm.",ifelse(!is.na(args$dropRep),paste0("drop",args$dropRep,"."),""),basename(args$HAFs))
cat("saving to ",filename,"\n")
save(df.glm,file=filename)


