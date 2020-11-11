
###############
## Utilities
##############

## redefine list function so you can return multiple variables from a function with one assignment call
#########################
list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

## chop a string and return a specified field
#########################
chop=function(x,splitchar,field) {sapply(strsplit(x,splitchar),"[",field)}

## get sizes of all objects in memory, sorted
get_obj_sizes=function(){print("sort( sapply(ls(),function(x){format(object.size(get(x)), units = 'Mb')}))")}

## apply a function to all pairs of columns in a matrix
fun.mat=function (x,myfun) {
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  
  ncy <- ncx <- ncol(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      x2 <- x[, i]
      y2 <- x[, j]
      r[i, j] <- myfun( x2, y2)
    }
  }
  r <- r + t(r) - diag(diag(r))
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
}


###############
## Transformations
###############

## calculate N-effective (Feder 2012), the variance of frequency estimates,
## approximated by calculating the effective number of observations at a given locus, conditional on read depth and number of chromosomes in the sample
n_eff=function(rd,pooledChroms){
  round((rd*pooledChroms-1)/(rd+pooledChroms))
}

### simulate panmixia given a vector of allele freq and matrix of read depths
create_panmixia=function(afmeans,rd){
  afpan=apply(rd,2,function(rd_i){
    rbinom(length(rd_i),rd_i,afmeans)/rd_i
  }
  )
}

## Recode allele freqs from 'ref-major-minor' to 'ref-alt'
recode_AF=function(df,afcols=NULL,coefcols=NULL){
  dropMe= (df$ref != df$major) & (df$ref != df$minor)
  df=df[!dropMe,]
  flipMe=df$ref == df$minor
  if(!is.null(afcols)){
    df[flipMe,afcols]=1-df[flipMe,afcols,with=F]
  }
  if(!is.null(coefcols)){
    df[flipMe,coefcols]=-1*df[flipMe,coefcols,with=F]
  }
  colnames(df)[grep("minor",colnames(df))]="alt"
  df$alt[flipMe]=df$major[flipMe]
  df=df[,-grep("major",colnames(df)),with=F]
  return(df)
}

## assign snps a folded freq bin based on alt allele freqs in founder set
assign_SNPs_freqBins=function(sites,snpFile,nFreqBins=10){
  #################
  ## sites is a dataframe listing the position of all SNPs, containing (at minimum) two columns: chrom,pos
  ## snpFile is the name of a text file containing the call (-1/0/0.5/1 corresponding to missing/ref/het/alt calls) in the founder set
  
  freqBins=do.call(rbind,lapply(chroms,function(chrom){
    snps=fread(gsub("CHROM",chrom,snpFile));snps[snps<0]=NA
    fb=rowSums(snps[,-1],na.rm=T)/rowSums(!is.na(snps[,-1]))
    fb[fb>.5]=1-fb[fb>.5]
    fb=cut(fb,seq(0,.5,by=.5/nfreqBins),labels=FALSE)*(.5/nfreqBins)
    return(data.frame(chrom=chrom,pos=snps[[1]],freq=fb))
  })) %>% merge(sites %>% dplyr::select(chrom,pos)) %>% arrange(chrom,pos) %>% mutate(freq=factor(freq))
}
############


###############
### Statistics
################

## get standard error
se <- function(x) sqrt(var(x)/length(x))

## calculate Fst between two (equi-length vectors of) allele freqs
Fst=function(p1, p2) {
  fhat <- p1/2 + p2/2  # avg freq across both pops ie (p1+p2)/2
  Htot <- 2*fhat*(1-fhat) #heterozygosity of the avg freq ie 2pq 
  Hwith <- p1*(1-p1) + p2*(1-p2) #avg heterozygosity of indivdual pop freqs  ie. (2pq + 2pq)/2
  fst=(Htot-Hwith)/Htot # how different are they? scaled by total
  fst[fhat==0 | fhat==1]=0 # set fst for fixed sites to zero
  return(fst)
}

## calculate Fst between all pairs of columns in an allele freqs matrix
Fst.mat=function (x) {
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  
  ncy <- ncx <- ncol(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq_len(ncx)) {
    for (j in seq_len(i)) {
      x2 <- x[, i]
      y2 <- x[, j]
      r[i, j] <- mean(Fst( x2, y2),na.rm=T)
    }
  }
  r <- r + t(r) - diag(diag(r))
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
}

## calculate Fst between more than two (equi-length vectors of) allele freqs
Fstmulti <- function(pAll) { #pAll should be a matrix with a separate column for each population
  fhat <- rowMeans(pAll)  # avg freq across pops 
  Htot <- 2*fhat*(1-fhat) #heterozygosity of the avg freq ie 2pq 
  Hwith <- rowMeans(2*pAll*(1-pAll)) #avg heterozygosity of indivdual pop freqs  
  (Htot-Hwith)/Htot # how different are they? scaled by total
}

## get zscore
zscore <- function(x,m=mean(x),s=sd(x)){
  z<- (x - mean(x)) / sd(x)
  return(z)
}

## calculate Heterozygosity from an allele frequency vector
hetz <- function(p) { 2*p*(1-p) }

## function to get drift-effective population size, following the method of Reed et al https://www.genetics.org/content/197/2/781.full
##--> based on expectation that the variance in AF over time is proportional to the frequency it started at and the number of gens since starting
get_drift_Ne=function(pInit,pEv,gen){
  
  flipMask=pInit>.5
  pInit[flipMask]=1-pInit[flipMask]
  pEv[flipMask]=1-pEv[flipMask]
  
  snpAFBinID=cut(pInit,breaks = seq(0,.5,.05),labels = FALSE)
  do.call(rbind,
          lapply(1:10,function(bin){
            mask=snpAFBinID==bin & pInit>.01 & pEv>.01 & pEv<.99
            nSites=sum(mask)
            if(nSites==0){
              return(data.frame(bin,initialMAF=NA,nSites,medianShift=NA,Ne=NA))
            } else{
              medianShift=median(abs(pEv[mask]-pInit[mask]))
              sqDev=(pEv[mask] - pInit[mask])^2/(pEv[mask]*(1-pEv[mask]))
              Ne=-gen/(2*log(1-1/nSites)*sum(sqDev)) 
              return(data.frame(bin,initialMAF=median(pInit[mask]),nSites,medianShift,Ne))
            }
          }))
}

### do PCA for sites on each chromosome and for genome
do_pca_perChrom=function(afmat,sites,centerMe=TRUE,scaleMe=TRUE,siteMask=NULL,chroms=NULL){
    if(is.null(chroms)){chroms=c("2L","2R","3L","3R","X","genome")}
    if(is.null(siteMask)){siteMask=rep(TRUE,nrow(sites))}
    list.pca=lapply(chroms,function(chrom){
      chrMask=if(chrom=="genome"){rep(TRUE,nrow(sites))} else{sites$chrom==chrom}
      return(prcomp(t(afmat[siteMask & chrMask,]),center = centerMe,scale=scaleMe))
    })
    names(list.pca)=chroms
    return(list.pca)
}

do_Fst_perChrom=function(afmat,sites,siteMask=NULL,chroms=NULL){
    if(is.null(chroms)){chroms=c("2L","2R","3L","3R","X","genome")}
    if(is.null(siteMask)){siteMask=rep(TRUE,nrow(sites))}
    list.fst=lapply(chroms,function(chrom){
      chrMask=if(chrom=="genome"){rep(TRUE,nrow(sites))} else{sites$chrom==chrom}
      return(Fst.mat(afmat[siteMask & chrMask,]))
    })
    names(list.fst)=chroms
    return(list.fst)
}


## function to fit a generalized linear model to allele frequencies at each site in siteIX (assumes afmat is already loaded)
## calculates the coefficient and pvalue for the contribution to the mdoel of each sample feature specified in model.vars
## can choose whether to fit model using 'glm' and quasibinomial error model, or 'glmer' with binomial error model, but the ability to specify some variables as random effects
## when using glmer can specify one or more sample features to be random effect variables (usually 'cage')
## specify the name of a sample feature in 'cmpAll' to convert that feature to a factor( if not already) and calculate coeffs/pvals for every pairwise comparison of factor levels
## specify the name of a sample feature in 'dont report' to omit coeffs/pvals from that feature from the returned results
fit_GLM=function(siteIX,sampIX,samps,model.vars,poolCt=100,glmType="glm",randEff=NULL,cmpAll=NULL,dontReport=NULL){
  df=samps[sampIX,colnames(samps)%in%model.vars,with=F];
  if(!is.null(cmpAll)){df[[cmpAll]]=factor(df[[cmpAll]])}
  
  formulaString=paste0(colnames(df),collapse=" + ")
  formulaString=paste0("cts ~ ",formulaString," + 1")
  
  if(glmType=="glmer"){
    sapply(randEff,function(rE){
      formulaString=gsub(rE,paste0('(1|',rE,')'),formulaString)
    })
  }
  cat("Model Forumla is: \n",formulaString)
  
  mclapply(siteIX,function(ix){
    
    Neff=((poolCt*2*rd[ix,sampIX])-1)/(poolCt*2+rd[ix,sampIX])
    cts=cbind(round(Neff*afmat[ix,sampIX]),round(Neff*(1-afmat[ix,sampIX])))
    
    model=switch(glmType,
                 glm=glm(as.formula(formulaString),family="quasibinomial",data=df),
                 glmer=glmer(as.formula(formulaString),family="binomial",nAGQ = 0,data=df)
    )
    if(!is.null(cmpAll)){
      model.multcomp=summary(eval(parse(text=paste0("glht(model, mcp(",cmpAll,"='Tukey'))")))) 
      cp=cbind(coefficients(model.multcomp),model.multcomp$test$pvalues)
      row.names(cp)=gsub(" - ","_",row.names(cp))
    } else{cp=summary(model)$coefficients[-1,c(1,4)]}
    colnames(cp)=c("coef","p")
    cp=cp[grep(dontReport,row.names(cp),invert = TRUE),,drop=F]
    results=c(cp[,1],cp[,2]);names(results)=c(paste0("coef.",row.names(cp)),paste0("p.",row.names(cp)))
    return(results)
  },mc.preschedule=TRUE)
  
}

set_up_sampData=function(sampIX,samps,model.vars,cmpAll){
  sampData=samps[sampIX,colnames(samps)%in%model.vars,drop=FALSE];
  sampData=sampData[,apply(sampData,2,function(x){length(unique(x))>1}),drop=F]
  if(!is.null(cmpAll)){sampData[[cmpAll]]=factor(sampData[[cmpAll]])}
  return(sampData)
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

fit_GLM_one = function(af.site,rd.site,sampData,formulaString,cmpAll=NULL,dontReport=NULL){
  Neff=((poolCt*2*rd.site)-1)/(poolCt*2+rd.site);
  cts=cbind(round(Neff*af.site),round(Neff*(1-af.site)))
  model=glm(as.formula(formulaString),family="quasibinomial",data=sampData) 
  cp=extract_coef_pval(model,cmpAll,dontReport)      
}

fit_GLM_all=function(siteIX,sampIX,samps,model.vars,poolCt=100,cmpAll=NULL,dontReport=NULL){
  
  sampData=set_up_sampData(sampIX,samps,model.vars,cmpAll)
  
  formulaString=paste0(colnames(sampData),collapse=" + ")
  formulaString=paste0("cts ~ ",formulaString," + 1")
  
  cat("Model Formula is: \n",formulaString,"\n")
  nSites=length(siteIX);cat("there are ",nSites," sites\n")
  do.call(rbind,mclapply(siteIX,function(ix){
    
    if(ix%%10000 == 0){cat("working on site ",ix,"\n")}
    
    cp=fit_GLM_one(afmat[ix,sampIX],rd$EC[rd$chrom==sites$chrom[ix]],sampData,formulaString,cmpAll,dontReport)
    results=c(cp[,1],cp[,2]);
    names(results)=c(paste0("coef.",row.names(cp)),paste0("p.",row.names(cp)))
    return(results)},mc.cores=ncores))
}

###################
## Annotations
##################
write_vcf=function(df,filename){
  colnames(df)[match(c("chrom","pos","ref","alt"),colnames(df))]=c("#CHROM","POS","REF","ALT")
  df$ID=1:nrow(df)
  df$QUAL=NA
  df$FILTER="PASS"
  df$INFO=NA
  cat("##fileformat=VCFv4.1\n",file = filename)
  suppressWarnings(write.table(df[,c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"),with=F],sep="\t",file=filename,row.names=F,quote=F,append = TRUE))
}

annotate_vcf=function(vcf,annotationScript,gffFeaturesList,gffFiles,mergeByPos=TRUE){
  system(paste0(annotationScript," ",vcf," ",gffFeaturesList,' "',gffFiles,'"'));
  ann=fread(gsub(".vcf$",".genes",vcf))
  if(mergeByPos){
    geneRows=ann$feature=="gene"
    ann.genenameagg=aggregate(ann$Name[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.genefullagg=aggregate(ann$fullname[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.flybaseagg=aggregate(ann$ID[geneRows],by=list(ann$chrom[geneRows],ann$pos[geneRows]),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann.featureagg=aggregate(paste(ann$Name,ann$feature,sep="|"),by=list(ann$chrom,ann$pos),FUN=function(x){paste(unique(x),collapse=";")})
    ann.entrezagg=aggregate(ann$EntrezGene,by=list(ann$chrom,ann$pos),FUN=function(x){paste(unique(x[!is.na(x)]),collapse=";")})
    ann=suppressWarnings(Reduce(function(x, y) merge(x, y, all=TRUE,by=c("Group.1","Group.2")), 
                                list(ann.genenameagg, ann.genefullagg,ann.flybaseagg,ann.featureagg,ann.entrezagg)
    ))
    colnames(ann)=c("chrom","pos","name","fullname","fbID","features","entrezID")
  }
  return(ann[order(ann$chrom,as.numeric(ann$pos)),])
}

add_snpEff=function(df,snpEffFile,fields="effect"){
  snpEff=fread(snpEffFile)
  df=cbind(df,snpEff[match(paste0(df$chrom,df$pos),paste0(snpEff$chrom,snpEff$pos)),fields,with=F])
  return(df)
}

add_goTerms=function(df,gotermData){
  load(gotermData)
  goMask=df.go$fbID%in%unlist(strsplit(df$fbID,";"))
  df.go=aggregate(df.go$goterm[goMask],by=list(df.go$fbID[goMask]),FUN=function(x){paste(x,collapse="|")})
  df$goterms=sapply(df$fbID,function(x){goMask=df.go$Group.1%in%strsplit(x,";");paste(df.go$x[goMask],collapse="|")})
  return(df)
}

mark_centro_telo_meres=function(df){  
  #from Comeron https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002905#s4
  #alternative: see Kauer et al 2003 http://www.genetics.org/content/165/3/1137 : loci with recombination rates <0.0001% recombination per kilobase after adjusting for zero recombination in males (i.e., multiplying by 0.67 for the X chromosome and by 0.5 for the third chromosome)
  chroms=c("2L","2R","3L","3R","X")
  starts=c(.5,5.2,0.7,9.7,2.3)*1000000
  ends=c(17.4,20.8,19.9,26.9,20.8)*1000000
  
  mask=do.call(c,mapply(function(chr,s,e){which(df$chrom==chr & (df$pos<s | df$pos>e))},chroms,starts,ends))
  return(seq_along(df$pos) %in% mask)
}


