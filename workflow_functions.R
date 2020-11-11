
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

library(data.table)
library(tidyr)



  
 ####################### 
 RunFullWorkflow=function(afFile,glmFile,snpFile,comparisons,cageSet,
                          fdrThreshs=c(.2,.05,.01),esThreshs=c(0,0.02,0.02),winSize=500,winShift=100,maxClusterBreak=100,
                          maxSNPPairDist=3000000,linkedClusterThresh=0.03,ncores=15){
##################
   ## load GLM data
    load(glmFile)
    sites=df.glm %>% dplyr::select(chrom,pos)
    
    ## find parallel sites at each timeseg
    FDR=get_glm_FDR(df.glm,comparisons)
    
    ## get af shifts in training cages and in test cage at each time segment
    load(afFile)
    afShifts=get_af_shifts(afmat,samps,cageSet,comparisons)
    
    pSig= (0 + (FDR<=fdrThreshs[1] & abs(afShifts)>=esThreshs[1])) + 
      (0 + (FDR<=fdrThreshs[2] & abs(afShifts)>=esThresh[2])) + 
      (0 + (FDR<=fdrThreshs[3] & abs(afShifts)>=esThresh[3]))
    
    df.sig=get_sig_sites(df.glm,comparisons,pSig,shuffle=FALSE)
    
    ### score windows
    df.wins=score_wins(df.sig ,sites,winSize,winShift) 
    df.sig.shuff=df.sig %>% group_by(timeseg) %>% mutate(ix=sample(1:nrow(sites),n())) %>%
      ungroup() %>% mutate(chrom=sites$chrom[ix],pos=sites$pos[ix])
    df.wins.shuff=score_wins(df.sig.shuff,sigLabels,sites,winSize,winShift) 
    df.winfdr=get_win_fdr(df.wins,df.wins.shuff)
    
    ### cluster windows and merge by linkage
    df.clust=cluster_wins(df.wins,df.winfdr,maxBreak=maxClusterBreak) %>% mutate(startPos=sites$pos[startSNP],endPos=sites$pos[endSNP]) 
    df.clust = df.sig %>% filter(sigLevel>1) %>% associate_snps_to_clusters(df.clust) %>% 
      find_snp_pairs(maxDist=maxSNPPairDist) %>% filter(pairType=="inter") %>%
      calc_Rsq_for_snp_pairs(ncores=ncores) %>% dplyr::select(-snp1.cl,-snp2.cl) %>%
      merge_linked_clusters(df.clust,df.sig,Rsq.thresh = linkedClusterThresh) 
    #save(df.clust,file=paste0("Rdata/GLMLOO.drop",dropCage,"_clusters.Rdata"))
    params=list(glmFile,comparisons,timesegs,cageSet,fdrThreshs,esThreshs,sigLabels,winSize=500,winShift=100,maxClusterBreak=100,maxSNPPairDist=3000000,linkedClusterThresh=0.03,ncores=15)
    results=list("sigSites"=df.sig,"wins"=df.wins,"clusters"=df.clust,"params"=params)
    return(results)
 }
 
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
 get_sig_sites=function(df.glm,comparisons,pSig,shuffle=FALSE){
   ########################
   if(shuffle){pIX=sample(1:nrow(df.glm))} else{pIX=1:nrow(df.glm)}
   do.call(rbind,lapply(comparisons,function(cc){
     cc.ix=match(cc,comparisons)
     df.glm %>% mutate(ix=1:nrow(df.glm)) %>% 
       rename(coef.div=paste0("coef.",cc),p.div=paste0("p.",cc)) %>% 
       select(ix,chrom,pos,coef.div,p.div) %>% 
       mutate(sigLevel=pSig[pIX,cc.ix],p.div=p.div[pIX],comparison=cc) %>% 
       filter(sigLevel>0,!is.na(p.div))
   })) %>% mutate(comparison=factor(comparison,comparisons)) 
 }
 ##########

 ########
 score_wins = function(df.sig,sites,winSize,winShift){
########################
    nSites=nrow(sites)
    comparisons=levels(df.sig$comparisons)
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
    comparisons=levels(df.wins$comparisons)
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
      
      
      
    },as.character(ct_pairs[,1]),as.character(ct_pairs[,2]),SIMPLIFY=FALSE)) %>%
      mutate(comparison=factor(comparison,comparisons))
  }
 ############
  
 ################
  associate_snps_to_clusters = function(df.snps,df.clust){
    do.call(rbind,lapply(1:nrow(df.clust),function(cl){
      df.snps %>% merge(df.clust[cl,] %>% dplyr::select(cl,comparison,chrom,startPos,endPos),by=c("comparison","chrom")) %>% 
        filter(pos>=C$startPos,pos<=C$endPos) 
    })) 
  }
 #############
  
 ################
 find_snp_pairs=function(df.snps,maxDist){
  ################
    do.call(rbind,lapply(unique(df.snps$chrom),function(chrom){
      do.call(rbind,lapply(unique(df.snps$comparison),function(cc){
        sigIX=which(df.snps$chrom==chrom & df.snps$comparison==cc)
        sigIX=expand.grid(sigIX,sigIX) %>% filter(Var1<Var2) 
        df.pairs=data.frame(snp1.pos=df.snps$pos[sigIX$Var1],snp2.pos=df.snps$pos[sigIX$Var2])
        if("cl" %in% colnames(df.snps)){
          df.pairs <- df.pairs %>% mutate(snp1.cl=df.snps$cl[sigIX$Var1],snp2.cl=df.snps$cl[sigIX$Var2]) %>%
            mutate(pairType=ifelse(snp1.cl==snp2.cl,"intra","inter"))
        }
        if("freq" %in% colnames(df.snps)){
          df.pairs <- df.pairs %>% mutate(snp1.freq=df.snps$freq[sigIX$Var1],snp2.freq=df.snps$freq[sigIX$Var2])
        }
        df.pairs <- df.pairs %>% mutate(snpDist=abs(snp1.pos-snp2.pos)) %>%
          filter(snpDist<maxDist) %>% mutate(chrom=chrom,timeseg=ts) 
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
        suppressWarnings(cor(t(snps[ii,-1]),t(snps[jj,-1]),use="p"))^2
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
        summarise(chrom=unique(chrom),comparison=unique(comparisons), nWins=sum(nWins),
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
  
  
  
#################
AssignSNPsFreqBins=function(sites,snpFile,nFreqBins=10){
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