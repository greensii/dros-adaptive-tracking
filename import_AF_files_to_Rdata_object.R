## make af data object
make_afData=function(thresh.altrd,thresh.minrd){
  if(check_files(afFiles)){download_files(afFiles)}
  if(check_files(repeatMaskerFiles)){download_files(repeatMaskerFiles)}
  
  samps=load_sampKey(sampleKeyFile) 
  afsamps=chop(basename(fread(system(paste0("ls ",afFiles[1]," | grep .af.meanrd$"),intern=T)[1])$V1),"_",1)
  sMask=match(samps$sampID,afsamps)
  
  cat("**Filtering sites from orchard_2017 dataset \n",file=siteFilterLog,append=FALSE)  
  #read af and filter
  df.af=foreach(chrom=chroms,.combine=rbind)%do%{
    rm(contammask,roundmask,foundermask, repmask,rdmask,refmask)
    
    af=fread(system(paste0("ls ",afFiles[1]," | grep ",chrom,".af$"),intern=T))
    sites=af[,1:5,with=F]
    rd=as.matrix(af[,seq(7,ncol(af),2),with=F])
    af=as.matrix(af[,seq(6,ncol(af),2),with=F])
    colnames(rd)=afsamps
    colnames(af)=afsamps
    cat(nrow(sites),",chrom",chrom,",chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter sites with potential human contamination from barcode swapping
    contam=fread(humanContamFilter[1])
    contammask=!(sites$pos%in%contam$pos[contam$chrom==chrom])
    sites=sites[contammask,];af=af[contammask,];rd=rd[contammask,]
    cat(sum(!contammask),",humanreads,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter out sites with non-overlapping allele freqs between round1 samps and round2 samps
    R1samps=afsamps%in%samps$sampID[samps$sequencingRd=="R1"]
    R1max=apply(af[,R1samps],1,max)
    R1min=apply(af[,R1samps],1,min)
    R2max=apply(af[,!R1samps],1,max)
    R2min=apply(af[,!R1samps],1,min)
    batchmask=!((R1min>R2max) | (R2min>R1max))
    sites=sites[batchmask,];af=af[batchmask,];rd=rd[batchmask,];
    cat(sum(!batchmask),",batcheffects,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## subset to relevant samps
    af=af[,sMask];rd=rd[,sMask];
    
    ## filter out sites not present in the baseline samples
    foundermask=apply(af[,samps$treatment=="Founder"],1,min)>thresh.founder
    sites=sites[foundermask,];af=af[foundermask,];rd=rd[foundermask,]
    cat(sum(!foundermask),",notstanding variation,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## filter repeat regions
    df.repeat<-fread(system(paste0("ls ",repeatMaskerFiles[1]," | grep ",chrom,".txt$"),intern=T))
    pos=sites$pos
    cc=cut(df.repeat$genoStart,c(0,pos),labels=FALSE)
    dd=cut(df.repeat$genoEnd,c(0,pos),labels=FALSE)
    cc[is.na(cc)]=max(cc,na.rm=TRUE)+1
    dd[is.na(dd)]=max(dd,na.rm=TRUE)+1
    pos.filter=eval(parse(text=paste0("c(",paste0(cc[dd>cc],":",dd[dd>cc]-1,collapse=","),")")))
    repmask=!(1:length(pos))%in%pos.filter
    sites=sites[repmask,];af=af[repmask,];rd=rd[repmask,]
    cat(sum(!repmask),",repeatregions,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    # filter out low + high rd sites
    ### since there are a few samps with low overall rd, only examine samps with mean overall rd>thresh.sampmeanrd
    sampsToExamine=colMeans(rd)>thresh.sampmeanrd
    rdmask=rowSums(rd[,sampsToExamine]<thresh.minrd)==0 & rowSums(rd>thresh.maxrd)==0 
    sites=sites[rdmask,];af=af[rdmask,];rd=rd[rdmask,]
    cat(sum(!rdmask),",lowreaddepth,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    ## keep only sites where ref allele is present
    refmask=(sites$ref==sites$major) | (sites$ref==sites$minor)
    sites=sites[refmask,];af=af[refmask,];rd=rd[refmask,]
    cat(sum(!refmask),",refallelenotfound,chrom",chrom,"\n",file=siteFilterLog,append=TRUE)
    
    cbind(sites,af,rd)
  }
  
  #save as Rdata
  sites=df.af[,1:5,with=F]
  afmat=as.matrix(df.af[,5+1:nrow(samps),with=F])
  rd=as.matrix(df.af[,5+nrow(samps)+1:nrow(samps),with=F])
  flipMe=(sites$ref != sites$major)
  afmat[flipMe,]=1-afmat[flipMe,]
  sites$minor[flipMe]=sites$major[flipMe]
  colnames(sites)[5]="alt";sites=sites[,-4,with=F]
  save(sites,samps,afmat,rd,file=afData)
  upload_files(afData)
}