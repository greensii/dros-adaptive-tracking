
## workflow steps:
#  get_glm_pvals
#  get_af_shifts
#  get_sig_sites
#  score_wins
#  get_win_fdr
#  cluster_wins
#  associate_snps_to_clusters
#  find_snp_pairs
#  calc_Rsq_for_snp_pairs
#  merge_linked_clusters

## additional functions:
#  get_loo_shifts
#  get_matched_pairs
#  get_af_trajectories

  
####################### 
RunFullWorkflow=function(afFile,glmFile,snpFile,
                        comparisons=NULL,cageSet=NULL,
                        fdrThreshs=c(.2,.05,.01),esThreshs=c(0,0.02,0.02),winSize=500,winShift=100,maxClusterBreak=100,
                        maxSNPPairDist=3000000,linkedClusterThresh=0.03,ncores=1){
##################
   ## load GLM data
   cat("loading GLM results..\n");flush.console();Sys.sleep(1)
    if(grepl(".rds$",tolower(glmFile))){df.glm=readRDS(glmFile)} else{
      load(glmFile)
    }
  sites=df.glm %>% dplyr::select(chrom,pos)
  
  ## select comparisons
  p_labs=gsub("^p.","",colnames(df.glm)[grepl("^p\\.",colnames(df.glm))])
  coef_labs=gsub("^coef.","",colnames(df.glm)[grepl("^coef\\.",colnames(df.glm))])
  if(!is.null(comparisons)){comparisons=intersect(intersect(comparisons,p_labs),coef_labs)}else{
    comparisons=intersect(p_labs,coef_labs)
  }
  cat("running analysis for comparisons:\n",paste0(paste0("+ ",comparisons),collapse="\n"));flush.console();Sys.sleep(1)
  
  ## load HAFs data
    cat("loading HAFs data..\n");flush.console();Sys.sleep(1)
    load(afFile)
    if(is.null(cageSet)){cageSet=unique(samps$cage)} else {
      cageSet=intersect(cageSet,unique(samps$cage))
    }
    cat("running analysis for cages:\n",paste0(paste0("+ ",cageSet),collapse="\n"));flush.console();Sys.sleep(1)
    
    ## get af shifts in training cages and in test cage at each time segment
    afShifts=get_af_shifts(afmat,samps,cageSet,comparisons)
    
    ## get FDR-corrected pvalues
    cat("finding parallel sites..\n");flush.console();Sys.sleep(1)
    FDR=get_glm_FDR(df.glm,comparisons)
    
    ## classify sites
    df.sig=get_sig_sites(df.glm,comparisons,FDR,afShifts,fdrThreshs,esThreshs)
    
    ### score windows
    cat("scoring windows..\n");flush.console();Sys.sleep(1)
    df.wins=score_wins(df.sig ,sites,winSize,winShift) 
    df.sig.shuff=df.sig %>% group_by(comparison) %>% mutate(ix=sample(1:nrow(sites),n())) %>%
      ungroup() %>% mutate(chrom=sites$chrom[ix],pos=sites$pos[ix])
    df.wins.shuff=score_wins(df.sig.shuff,sites,winSize,winShift) 
    df.winfdr=get_win_fdr(df.wins,df.wins.shuff)
    
    ### cluster windows and merge by linkage
    cat("finding initial clusters..\n");flush.console();Sys.sleep(1)
    df.clust=cluster_wins(df.wins,df.winfdr,maxBreak=maxClusterBreak) %>% mutate(startPos=sites$pos[startSNP],endPos=sites$pos[endSNP]) 
    
    cat("merging linked clusters..\n");flush.console();Sys.sleep(1)
    df.clust = df.sig %>% filter(sigLevel>1) %>% associate_snps_to_clusters(df.clust) %>% 
      find_snp_pairs(maxDist=maxSNPPairDist) %>% filter(pairType=="inter") %>%
      calc_Rsq_for_snp_pairs(ncores=ncores,snpFile) %>% dplyr::select(-snp1.cl,-snp2.cl) %>%
      merge_linked_clusters(df.clust,df.sig,Rsq.thresh = linkedClusterThresh) 
    
    cat("extracting trajectories for top cluster sites..\n");flush.console();Sys.sleep(1)
    df.traj=get_af_trajectories(chrompos=df.clust %>% select(chrom,bestSNP.pos),HAFsFile,baselineFile)
    
    cat("finished.\n");flush.console();Sys.sleep(1)
    cat(nrow(df.sig),"parallel sites in",nrow(df.clust),"clusters\n");flush.console();Sys.sleep(1)
    
    params=list(glmFile,comparisons,cageSet,fdrThreshs,esThreshs,
                winSize,winShift,maxClusterBreak,maxSNPPairDist,linkedClusterThresh,ncores)
    results=list("sigSites"=df.sig,"wins"=df.wins,"clusters"=df.clust,"clMarkerAF"=df.traj,"params"=params)
    return(results)
 }
########################


#############
 get_glm_FDR=function(df.glm,comparisons){
################
   pvals=as.matrix(df.glm %>% dplyr::select(paste0("p.",comparisons)))
   pvals[pvals==0]=1;
   FDR=matrix(p.adjust(pvals,method="BH"),nrow(pvals),ncol(pvals))
   colnames(FDR)=paste0("p.",comparisons)
   return(FDR)
 }
 #################
 
 #################
 get_af_shifts=function(afmat,samps,cage_set=NULL,comparisons){
   ################
   df.shifts=do.call(cbind,lapply(comparisons,function(cc){
     tt=rev(chop(cc,"_",1:2));
     if(is.null(cage_set)){cage_set=unique(samps$cage)}
     cageMask=samps$cage %in% cage_set
     t1Mask=cageMask & samps$tpt==tt[1]; 
     if(sum(t1Mask)>1){af1=rowMeans(afmat[,t1Mask])} else{af1=afmat[,t1Mask]}
     t2Mask=cageMask & samps$tpt==tt[2]
     if(sum(t2Mask)>1){af2=rowMeans(afmat[,t2Mask])} else{af2=afmat[,t2Mask]}
     ss=data.frame(dAF=af2-af1) %>% 
       rename_all(.funs=function(x){paste0("dAF.",cc)})  
     return(ss)
   }))
   return(df.shifts)
 }
 ################
 
 #############
 get_sig_sites=function(df.glm,comparisons,FDR,afShifts,fdrThreshs,esThreshs){
   ########################
  
   pSig= Reduce("+",lapply(1:length(fdrThreshs),function(ii){
     (0 + (FDR<=fdrThreshs[ii] & abs(afShifts)>=esThreshs[ii]))
   }))
   
   do.call(rbind,lapply(comparisons,function(cc){
     cc.ix=match(cc,comparisons)
     df.glm %>% mutate(ix=1:nrow(df.glm)) %>% 
       rename(coef.div=paste0("coef.",cc),p.div=paste0("p.",cc)) %>% 
       select(ix,chrom,pos,coef.div,p.div) %>% 
       mutate(sigLevel=pSig[,cc.ix],FDR=FDR[,cc.ix],afShift=afShifts[,cc.ix],comparison=cc) %>% 
       filter(sigLevel>0,!is.na(p.div))
   })) %>% mutate(comparison=factor(comparison,comparisons)) 
 }
 ##########

 ########
 score_wins = function(df.sig,sites,winSize,winShift){
########################
    nSites=nrow(sites)
    chroms=unique(df.sig$chrom)
    comparisons=levels(df.sig$comparison)
    
    ## get average sigfreq per base, with sig weighted by sig level
    sigCounts=df.sig %>% group_by(sigLevel) %>% summarise(count=n())
    #print(sigCounts); flush.console()
    
    cat(winSize,"-SNP windows with ",winShift,"-SNP shift\n")
    
    ## move along windows and get scores
    winStarts=seq(1,nSites,by=winShift);
    
    df.null=data.frame(expand.grid(chroms,comparisons),stringsAsFactors = FALSE) %>% 
    setNames(c("chrom","comparison")) %>% mutate_all(.funs = as.character) %>% mutate(sigLevel=0,p.div=1,pos=0)
    df.win=do.call(rbind,mclapply(winStarts,function(ww){
      winChrom=sites$chrom[ww];
      #if(ww%%1000==0){cat(winChrom,", snp ",ww,"\n"); flush.console()}
      if((ww+winSize)<nSites & (winChrom==sites$chrom[ww+winSize])){
        df.sig %>%  filter(ix>=ww,ix<ww+winSize) %>% 
          select(chrom,pos,comparison,sigLevel,p.div) %>% 
          rbind(df.null %>% dplyr::filter(chrom==winChrom)) %>% 
          group_by(chrom,comparison) %>% summarise(winscore=sum(sigLevel),
                                                nSig.1=sum(sigLevel==1),
                                                nSig.2=sum(sigLevel==2),
                                                nSig.3=sum(sigLevel==3),
                                                bestP=ifelse((nSig.2+nSig.3)>0,min(p.div[sigLevel>1]),1),
                                                bestSNP=ifelse((nSig.2+nSig.3)>0,pos[bestP==p.div],0)
          ) %>% 
          mutate(winstartSNP=ww,winstopSNP=ww+winSize,winSize=winSize) %>% ungroup()
      } else(return(NULL))
    },mc.cores=12)) %>% mutate(comparison=factor(comparison,comparisons))
    
    return(df.win)
  }
 ########

 #############
 get_win_fdr=function(df.wins,df.wins.shuff){
 ############
    do.call(rbind,lapply(levels(df.wins$comparison),function(cc){
      ss=df.wins.shuff %>% filter(comparison==cc) %>% 
        group_by(winscore) %>% summarise(count=n()) %>% 
        ungroup() %>% arrange(desc(winscore))
      
      fdr=0;ii=0
      while(fdr<.05){
        ii=ii+1
        score.shuff=ss$winscore[ii];cat(score.shuff,"..");
        fdr=sum(ss$count[1:ii])/sum(df.wins$comparison==cc & df.wins$winscore>=score.shuff)
      };cat("\n")
      return(data.frame(comparison=cc,ii=ii,score=score.shuff,fdr,
                        baserate=(df.wins.shuff %>% filter(comparison==cc) %>% pull(winscore) %>% mean)
                        )
          )
    }))
 }
 #############
  
 ################
 cluster_wins = function(df.wins,df.winfdr,maxBreak){
 ########################
    
    ##filter for windows with a score above sigRate
    df.wins=df.wins %>%  arrange(chrom,comparison,winstartSNP,winSize) 
    chroms=unique(df.wins$chrom)
    comparisons=levels(df.wins$comparison)
    cc_pairs=expand.grid(chroms,comparisons)
    
    ## for each chrom/comparison, go through windows one by one - 
    ## if start of one window is before end of current cluster (plus maxBreak distance), add it to the cluster and extend cluster end
    df.clust = do.call(rbind,mcmapply(function(chr,cc){
      baserate=df.winfdr %>% filter(comparison==cc) %>% pull(baserate)
      sigrate=df.winfdr %>% filter(comparison==cc) %>% pull(score)
      cat(chr,cc,baserate,sigrate,"\n")
      ww= df.wins %>% filter(chr==chrom & cc==comparison & winscore>=baserate) 
      if(nrow(ww)==0){return(NULL);cat("no wins\n")}
      
      ww <- ww %>% mutate(cl=NA)
      
      cl=0;cStart=0;cEnd=0
      for(ii in 1:nrow(ww)){ 
        if(ww$winstartSNP[ii]<=(cEnd+maxBreak)){
          ww$cl[ii]=cl;cEnd=max(cEnd,ww$winstopSNP[ii])
        } else{ 
          cl=cl+1;ww$cl[ii]=cl;cStart=ww$winstartSNP[ii];cEnd=ww$winstopSNP[ii]
        }
      }
      # cat(nrow(ww),"wins\n")
      ww <- ww %>% group_by(cl) %>% 
        summarise(chrom=unique(chrom),comparison=unique(comparison), nWins=n(),
                  startSNP=min(winstartSNP),endSNP=max(winstopSNP),
                  winscore.mean=mean(winscore),winscore.max=max(winscore),
                  nSig.1=sum(nSig.1),nSig.2=sum(nSig.2),nSig.3=sum(nSig.3),
                  bestSNP.p=min(bestP),bestSNP.pos=bestSNP[bestP==bestSNP.p][1]
        ) %>% ungroup() %>% filter((nSig.2+nSig.3)>0,winscore.max>=sigrate) 
      
      if(nrow(ww)>0){
        ww <- ww %>% mutate(cl=paste0("c",1:length(cl),".",chrom,".",comparison))
        # cat(nrow(ww),"clusters\n")
        return(ww)
      } else{
        #  cat(nrow(ww),"clusters\n");
        return(NULL)
      }
      
      
      
    },as.character(cc_pairs[,1]),as.character(cc_pairs[,2]),SIMPLIFY=FALSE)) %>%
      mutate(comparison=factor(comparison,comparisons))
  }
 ############
  
 ################
  associate_snps_to_clusters = function(df.snps,df.clust){
    do.call(rbind,lapply(1:nrow(df.clust),function(cl){
      df.snps %>% merge(df.clust[cl,] %>% dplyr::select(cl,comparison,chrom,startPos,endPos),by=c("comparison","chrom")) %>% 
        filter(pos>=startPos,pos<=endPos) 
    })) 
  }
 #############
  
 ################
 find_snp_pairs=function(df.snps,maxDist,cross_ts=FALSE){
  ################
    do.call(rbind,lapply(unique(df.snps$chrom),function(chrom){
      do.call(rbind,lapply(unique(df.snps$comparison),function(cc){
        sigIX1=which(df.snps$chrom==chrom & df.snps$comparison==cc)
        if(cross_ts){
          sigIX2=which(df.snps$chrom==chrom & df.snps$comparison!=cc)
        } else{
          sigIX2=which(df.snps$chrom==chrom & df.snps$comparison==cc)
        }
        sigIX=expand.grid(sigIX1,sigIX2) %>% filter(Var1<Var2) 
        df.pairs=data.frame(snp1.pos=df.snps$pos[sigIX$Var1],snp2.pos=df.snps$pos[sigIX$Var2])
        if("cl" %in% colnames(df.snps)){
          df.pairs <- df.pairs %>% mutate(snp1.cl=df.snps$cl[sigIX$Var1],snp2.cl=df.snps$cl[sigIX$Var2]) %>%
            mutate(pairType=ifelse(snp1.cl==snp2.cl,"intra","inter"))
        }
        if("freq" %in% colnames(df.snps)){
          df.pairs <- df.pairs %>% mutate(snp1.freq=df.snps$freq[sigIX$Var1],snp2.freq=df.snps$freq[sigIX$Var2])
        }
        df.pairs <- df.pairs %>% mutate(snpDist=abs(snp1.pos-snp2.pos)) %>%
          filter(snpDist<maxDist) %>% mutate(chrom=chrom,comparison=cc) 
        return(df.pairs)
        
      }))
    })) 
  }
 ################
  
 ################# 
  calc_Rsq_for_snp_pairs=function(snppairs,ncores,snpFile){
 ###############
    ## snppairs must have columns: chrom, snp1.pos, snp2.pos
    snppairs$Rsq=NA
    for(mychrom in sort(unique(snppairs$chrom))){
      cat(mychrom,"...")
      snps=fread(gsub("CHROM",mychrom,snpFile));snps[snps<0]=NA
      
      mask=snppairs$chrom==mychrom
      ps=snppairs[mask,] %>%
        mutate(
          snp1.ix=match(snp1.pos,snps[[1]]),
          snp2.ix=match(snp2.pos,snps[[1]])
        )
      cat(nrow(ps),"pairs...")
      snppairs$Rsq[mask]=mcmapply(function(ii,jj){
        Rsq=NA
        try({Rsq=suppressWarnings(cor(t(snps[ii,-1]),t(snps[jj,-1]),use="p"))^2})
        return(Rsq)
      },ps$snp1.ix,ps$snp2.ix,mc.cores=ncores)
      
      cat("done.\n")
    }
    return(snppairs)
  }
  ##################
  
 ##################
  merge_linked_clusters=function(snppairs,df.clust,df.sig,Rsq.thresh){
 ##################
    linkage_remaining=TRUE
    linked_clusters=NULL
    
    while(linkage_remaining){
      
      df.clust <- df.clust %>% mutate(clnew=cl)
      if(!is.null(linked_clusters)){
        for(ii in 1:nrow(linked_clusters)){
          c1=linked_clusters$snp1.cl[ii];c2=linked_clusters$snp2.cl[ii]
          if((c1 %in% df.clust$clnew) & (c2 %in% df.clust$clnew)){
            df.clust$clnew[match(c(c1,c2),df.clust$clnew)]=paste0(c1,"-",c2)
          }}
      }
      df.clust <- df.clust %>% group_by(clnew) %>%
        summarise(chrom=unique(chrom),comparison=unique(comparison), nWins=sum(nWins),
                  startSNP=min(startSNP),endSNP=max(endSNP),
                  startPos=min(startPos),endPos=max(endPos),
                  winscore.mean=mean(winscore.mean),winscore.max=max(winscore.max),
                  nSig.1=sum(nSig.1),nSig.2=sum(nSig.2),nSig.3=sum(nSig.3),
                  bestSNP.pos=bestSNP.pos[bestSNP.p==min(bestSNP.p)][1],bestSNP.p=min(bestSNP.p)
        ) %>% ungroup() %>% group_by(chrom,comparison) %>% arrange(startPos) %>%
        mutate(cl=paste0("c",1:length(clnew),".",chrom,".",comparison)) %>%
        ungroup()
      cat(nrow(df.clust)," clusters\n")
      
      
      df.sig.cl=associate_snps_to_clusters(df.sig %>% filter(sigLevel>1),df.clust) 
      cat(nrow(df.sig.cl),"/",nrow(df.sig %>% filter(sigLevel>1))," sig sites in clusters\n")
      
      
      linked_clusters=snppairs %>% 
        merge(df.sig.cl %>% dplyr::select(chrom,pos,comparison,cl) %>% rename(snp1.pos=pos,snp1.cl=cl),
              by=c("chrom","snp1.pos","comparison")) %>%
        merge(df.sig.cl %>% dplyr::select(chrom,pos,comparison,cl) %>% rename(snp2.pos=pos,snp2.cl=cl),
              by=c("chrom","snp2.pos","comparison")) %>%
        filter(!is.na(snp1.cl),!is.na(snp2.cl)) %>% mutate(cnum1=as.numeric(gsub("c","",chop(snp1.cl,"\\.",1))),
                                                           cnum2=as.numeric(gsub("c","",chop(snp2.cl,"\\.",1)))) %>%
        filter(cnum2==(cnum1+1))%>% 
        group_by(snp1.cl,snp2.cl) %>% summarise(ct=n(),Rsq.mean=mean(Rsq),Rsq.med=median(Rsq),Rsq.sd=sd(Rsq)) %>%
        filter(Rsq.mean>Rsq.thresh) %>% arrange(desc(Rsq.mean))
      
      cat(nrow(linked_clusters),"cluster pairs with high linkage\n")
      if(nrow(linked_clusters)==0){linkage_remaining=FALSE}
    }
    return(df.clust %>% dplyr::select(-clnew))
  }
 #################
  
  ##################
  ## ADDITIONS TO WORKFLOW
  ##################
  
 ########################
 get_af_trajectories=function(chrompos,HAFsFile,baselineFile){
########################
    list[sites,samps,afmat]=load_afs_with_baseline(HAFsFile,baselineFile)
    siteIX=sites %>% mutate(ix=1:nrow(sites)) %>% merge(chrompos,by.x=c("chrom","pos"),by.y=colnames(chrompos)) %>% pull(ix)
    
    samps <- samps %>% mutate(tpt=ifelse(is.na(tpt),0,tpt))
    
    df=do.call(rbind,lapply(siteIX,function(ix){
      samps %>% mutate(af=afmat[ix,],ix,chrom=sites$chrom[ix],pos=sites$pos[ix])
    }))
    
    return(df)
  }
 ######################
  
 ########################
 get_loo_shifts=function(results.loo,afFile,snpFile,comparisons){
  ########################
    load(afFile)
    freqBins=assign_SNPs_freqBins(sites,snpFile)
    shiftGroups=c("significant parallel shifts","significant anti-parallel shifts","shifts not significant")
    
    get_pval=function(set1,set2,sigMetric,paired){
      if(sigMetric=="ttest"){return(t.test(set1,set2,paired=paired)$p.value)}
      if(sigMetric=="ranksum"){return(wilcox.test(set1,set2,paired=paired)$p.value)}
    }
    
    shifts.loo=do.call(rbind,lapply(1:length(results.loo),function(dropCage){
      
      ## get af shifts for training cages and left-out cage
      shifts.train=get_af_shifts(afmat,samps,cage_set=cages[-dropCage],comparisons)
      shifts.test=get_af_shifts(afmat,samps,cage_set=dropCage,comparisons)
      
      ## find identified sig sites [ix,comparison] 
      sigIX=results.loo[[dropCage]][["sigSites"]] %>% filter(sigLevel>1) %>% 
        dplyr::select(ix,comparison) %>% mutate(comparison=as.numeric(comparison)) %>% as.matrix
      
      ## label list of site freqs as sig or not
      shifts=freqBins %>% mutate(ix=1:nrow(freqBins),sig=ifelse(ix%in% sigIX[,1],"sig","matched")) %>%
        ## for each chrom and freqbin, give a random index to each sig site and each non-sig site
        group_by(chrom,freq,sig) %>% mutate(rr=sample(rank(pos))) %>%
        ungroup() %>% dplyr::select(-pos) %>%
        ## spread list so sites with the same rr index are lined up
        spread(key="sig",val="ix") %>% 
        ## get rid of sites without a matching sig site
        filter(!is.na(sig)) %>% 
        ## add back info for the sig site
        merge(sigIX,by.x="sig",by.y="ix",all=T) %>% 
        ## turn back into a list, sig sites first, then matched
        gather(key="site",val="ix",sig,matched) %>%
        mutate(rr=paste0(rr,".",chrom,".",freq)) %>%
        mutate(trainDir=c(-1,1)[1+as.numeric(shifts.train[cbind(ix,comparison)]>0)]) 
      
      shifts<-cbind(shifts %>% dplyr::select(chrom,comparison,site,rr),
                    shifts.test[shifts$ix,] %>% mutate_all(function(x){x*shifts$trainDir}) 
                    )  %>% gather(key="c.meas",val="shift",-chrom,-comparison,-site,-rr) %>% 
        mutate(c.meas=match(chop(c.meas,"\\.",2),comparisons)) %>%
        mutate(dropCage=dropCage)   
  }))
    
    medians.loo = shifts.loo %>% arrange(rr) %>% group_by(comparison,dropCage,c.meas) %>% 
      mutate(pval.par.paired=get_pval(shift[site=="matched"],shift[site=="sig"],'ttest',TRUE),
             pval.nonpar.paired=get_pval(shift[site=="matched"],shift[site=="sig"],'ranksum',TRUE),
             pval.par.npaired=get_pval(shift[site=="matched"],shift[site=="sig"],'ttest',FALSE),
             pval.nonpar.npaired=get_pval(shift[site=="matched"],shift[site=="sig"],'ranksum',FALSE)) %>% 
      group_by(comparison,site,c.meas,dropCage) %>% 
      summarise(nSites=n(),nPos=sum(shift>0),shift=median(shift),
                pval.par.paired=unique(pval.par.paired),pval.nonpar.paired=unique(pval.nonpar.paired),
                pval.par.npaired=unique(pval.par.npaired),pval.nonpar.npaired=unique(pval.nonpar.npaired)) %>%
      ungroup() %>%
      mutate(shiftGroup=ifelse(p.adjust(pval.par.npaired,method = "BH")<.05,ifelse(shift>0,shiftGroups[1],shiftGroups[2]),shiftGroups[3])) %>%
      mutate(
        tpt.ID=factor(paste0("SNPs ID'd in 9 cages\n",comparisons[comparison]),
                      paste0("SNPs ID'd in 9 cages\n",comparisons)),
        tpt.measured=factor(paste0("shifts\n",comparisons[c.meas]),
                            paste0("shifts\n",comparisons))
      )
    return(list(shifts.loo,medians.loo))
  }
 ########################
 
 ########################
 get_matched_pairs=function(clusterPairs,snpFile,freqBins,distBinSize=25,maxPairsPerCluster=NULL,nCores=1){
   ####################
   #load("Rdata/glm.Ecages.Rdata");
   #sites=df.glm %>% dplyr::select(chrom,pos); rm(df.glm)
   #freqBins=assign_SNPs_freqBins(sites,snpFile) %>% mutate(ix=1:nrow(freqBins))
   freqBins=freqBins %>% mutate(ix=1:nrow(freqBins))
   clusterPairs <- clusterPairs %>% merge(freqBins,by.x=c("chrom","snp1.pos"),by.y=c("chrom","pos")) %>% rename(snp1.freq=freq,snp1.ix=ix) %>%
     merge(freqBins,by.x=c("chrom","snp2.pos"),by.y=c("chrom","pos")) %>% rename(snp2.freq=freq,snp2.ix=ix) 
   
   if(!is.null(maxPairsPerCluster)){
     clusterPairs <- clusterPairs %>% group_by(snp1.cl,snp2.cl) %>% mutate(rr=sample(order(snpDist))) %>% filter(rr<maxPairsPerCluster) %>% ungroup()
   }
   # for each pair of SNPs in clusterPairs (either pairs in the same cluster or pairs in different clusters)
   freqBins$sig=(freqBins$ix %in% clusterPairs$snp1.ix) | (freqBins$ix %in% clusterPairs$snp2.ix)
   
   system.time({candidatePairs=do.call(rbind,mclapply(1:nrow(clusterPairs),function(pp){
     #        cat(pp,"..");flush.console();Sys.sleep(1)
     #candidatePairs=do.call(rbind,mclapply(1:500,function(pp){
     ## get indexes of SNPs in frequency bins that match SNP1 and SNP2 but are not in clusters
     ix1=which(freqBins$chrom==clusterPairs$chrom[pp] & freqBins$freq==clusterPairs$snp1.freq[pp] & !freqBins$sig);
     ix2=which(freqBins$chrom==clusterPairs$chrom[pp] & freqBins$freq==clusterPairs$snp2.freq[pp] & !freqBins$sig);
     ## get the target distance 
     dB=abs(clusterPairs$snp1.pos[pp]-clusterPairs$snp2.pos[pp]);
     
     ## try to find a matched reference pair
     snp2.ix=NULL;
     while(length(snp2.ix)==0){
       ## pick a random SNP with freq matched to SNP1
       snp1.ix=sample(ix1,1)
       pos.1=freqBins$pos[snp1.ix];
       ## get the distance from that SNP to all SNPs with freq matched to SNP2
       dist.2=abs(freqBins$pos[ix2]-pos.1)
       ## pick pair with closest dist to target distance
       best=min(abs(dist.2-dB))
       if(best<=distBinSize){
         snp2.ix=ix2[abs(dist.2-dB)==best][1]
       }
       
     };
     
     ## return dataframe with found SNP sets
     data.frame(pp,snp1.ix,snp2.ix,pos.1) %>% mutate(pos.2=freqBins$pos[snp2.ix])
   },mc.cores = nCores))})
   return(candidatePairs)
 }
 #########################
  