## plotting_functions

#####################
## helper functions
####################
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 5, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

gg_color_hue=function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

plot_site=function(af,samps){
  samps$af=af
  ggplot(samps,aes(x=tpt,y=af,color=cage)) + geom_point() + facet_wrap(facets=vars(treatment),scales = "free_x") + 
    theme_minimal()
}

plot_MDS=function(df.fst,samps,treatments,tpts,drawArrows=FALSE,mytitle){
  sampMask=samps$treatment%in%treatments & samps$tpt%in%tpts
  fit<-cmdscale(as.matrix(df.fst),k=2,eig=TRUE)
  
  df=cbind(as.data.frame(fit$points[sampMask,]),samps[sampMask,])
  df$group=paste0(df$treatment,df$cage)
  
  df2=df[order(df$treatment,df$cage,df$tpt),]
  df2$V1.end=c(df2$V1[-1],0)
  df2$V2.end=c(df2$V2[-1],0)
  df2=df2[df2$tpt<4,]
  df2[df2$treatment=="FOUND",c("V1.end","V2.end")]=df2[df2$treatment=="FOUND",c("V1","V2")]
  df2$shift.y=df2$V2.end-df2$V2
  df2$shift.x=df2$V1.end-df2$V1
  
  ### color by treatment,alpha by timept, lines within cage
  ggMDS<-ggplot(df,aes(V1,V2,label=sampID,color=treatment)) + 
    labs(x="MDS1", y="MDS2", color="Treatment", alpha="Timepoint") + 
    geom_point(size=2,aes(alpha=tpt)) +
    ggtitle(mytitle) +
    scale_alpha(range = c(0.4, 1),limits=c(1,4),breaks=sort(unique(df$tpt)),labels=sort(unique(df$tpt))) +
    scale_color_manual(labels=c('Founder' ="Founder","E"="E","F"="N","ED"="ED","B"="B","S"="S"),
                       values = c('Founder' = "black", 
                                  'E' = gg_color_hue(3)[1],
                                  'F' = gg_color_hue(3)[2],
                                  'ED' = gg_color_hue(3)[3],
                                  'B' = gg_color_hue(5)[2],
                                  'S' = gg_color_hue(5)[5])) + 
    theme_minimal() +  
    theme(plot.title = element_text(size=12, face="bold"),axis.text = element_blank())
  
  if(drawArrows){
    ggMDS<- ggMDS + geom_segment(data=df2,
                                 aes(x=V1,y=V2,xend=V1.end,yend=V2.end,col=treatment,alpha=tpt+1,linetype="Cage"),
                                 arrow=arrow(type = "closed",length = unit(0.3,"cm")),
                                 show.legend = FALSE) +
      scale_linetype_manual(values=c("Cage" = 4))
  }
  
  return(ggMDS)
}

make_window_GLM_plot=function(df,scoreCol,colorCol,myalpha=.7,myshape=1){
  ### plots scores for window-based enrichment as a point in the center of the window
  df$score=df[[scoreCol]]
  df$color=df[[colorCol]]
  gg<-ggplot(df,aes((winStart+winStop)/2000000,score)) + 
    geom_point(aes(color=color),alpha=myalpha,shape=myshape) + 
    theme_minimal() + facet_wrap(~chrom,nrow = 1,scales="free_x") + 
    labs(x="genomic window (Mb)",y=expression("-log"[10] * "(p"["enrichment"] *")"))
  gg
}

###################
## main functions
##################

##################
make_loo_plot=function(medians.loo){
##################
  medians.loo %>% { ggplot(data=filter(.,site=="sig")) + 
      geom_jitter(aes(x=as.numeric(tpt.measured) -.2,y=shift,color=shiftGroup),width=.05,height=0,size=2,shape=19,alpha=.7) +
      geom_jitter(data=filter(.,site=="matched"),aes(x=as.numeric(tpt.measured) +.2,y=shift),color="darkslategray",width=.05,height=0,size=4,shape="x",alpha=.7)
  } + facet_wrap(~ tpt.ID,scales = "free_x",nrow=1) + 
    labs(x="",y="median phased shift in 10th cage",color="",shape="") +
    scale_color_manual(values=c("gray","lightgreen","coral")) + 
    scale_x_continuous(breaks=1:5,labels=levels(medians.loo$tpt.measured),expand = expand_scale(add = .4)) + 
    geom_hline(yintercept = 0) + #
    theme_minimal() + theme(panel.background  = element_rect(fill = 'transparent',color="darkslategray"),
                            panel.grid.major.x = element_blank(),
                            legend.position = "top",
                            strip.text = element_text(size=10),
                            axis.text.x=element_text(size=7)) + 
    lims(y=c(-.075,.075))
}


##################
make_manhattan=function(df.sig,df.clust=NULL,comparisonLabels=NULL,sigLevelLabels=NULL,sigColors=NULL,
                        show_inversions=TRUE,show_cluster_bars=FALSE,show_cluster_lines=FALSE,show_cluster_points=TRUE,show_cluster_nums=FALSE,
                        cluster_line_color="gray",ylims=c(-3,11)){
########################
  df.sig <- df.sig %>% mutate(sigLabel=factor(sigLevelLabels[sigLevel],sigLevelLabels),compLabel=factor(comparisonLabels[comparison],comparisonLabels))
  df.clust <-df.clust %>% group_by(comparison,chrom) %>% mutate(ct=n()) %>% mutate(clnum=rank(startPos)) %>% ungroup() %>%
    mutate(clfill=factor(clnum%%2)) %>% gather(key=posType,val=pos,startPos,endPos,bestSNP.pos) %>% 
    mutate(compLabel=factor(comparisonLabels[comparison],comparisonLabels)) 
  
  ggbase <- df.sig %>% ggplot(aes(x=pos/1000000))+ theme_minimal() 
  
  
  gg<- ggbase + lims(y=ylims)  +
    facet_grid(compLabel ~ chrom,scales="free_x", switch="y") + 
    theme(legend.position = "top",strip.text = element_text(size=16),strip.text.y = element_text(size=12)) +
    labs(x="Position (Mb)",y="-log10(parallelism p-value)")  + theme(strip.text.y = element_text(size=10))
  
  
  if(show_cluster_bars){
    gg<-gg+geom_ribbon(data =  df.clust, aes(group=cl,fill=clfill),ymin=0,ymax=Inf,na.rm=TRUE,alpha=.5) + 
      # scale_fill_gradientn(colours=brewer.pal(2,"Spectral"))
      scale_fill_viridis_d(option="cvidis",begin=.4,end = .7)  +
      guides(fill=FALSE) 
  }
  
  gg <- gg + suppressWarnings(geom_point(aes(y=-log10(p.div),color=sigLabel),alpha=.8,size=1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1,size=3))) +
    scale_color_manual(values=sigColors) + labs(color="")
  if(show_cluster_nums){
    gg<-gg+geom_text_repel(data = df.clust %>% filter(posType=="bestSNP.pos") ,
                           aes(y=0,label=clnum),size=3,alpha=1) 
  }
  
  if(show_cluster_lines){
    gg<-gg+
      geom_ribbon(data = df.clust ,
                  aes(fill=af0,group=cl),ymin=10,ymax=11,alpha=.49) 
    
  }
  
  if(show_cluster_points){
    gg<-gg+
      geom_point(data = df.clust %>% filter(posType=="bestSNP.pos"),
                 shape=19,y=0,color="darkslategray",size=.4) 
    
    
  }
  
  
  
  if(show_inversions){
    inv.breaks=inv.breaks %>% group_by(chrom) %>% mutate(rr=rank(start)) %>% ungroup() %>% 
      mutate(compLabel=factor(comparisonLabels[length(comparisonLabels)],comparisonLabels)) 
inv.lines= inv.breaks %>% ungroup() %>% gather(key=posType,val=pos,start,stop)
inv.labs=inv.breaks%>%group_by(chrom) %>% mutate(rr=rank(start),timeseg=ts_factor(timeseg.last)) %>% ungroup() %>% gather(key=winSide,val=pos,start,stop)%>%group_by(chrom,rr,name,timeseg)%>% summarise(pos=mean(pos))%>%ungroup()

gginv<- ggbase + suppressWarnings(geom_point(data=df.sig %>% group_by(chrom) %>% filter(pos==min(pos) | pos==max(pos)),
                                             aes(y=0),color="transparent")) +
  geom_line(data=inv.lines,aes(group=name,y=-rr),alpha=1,color="black") + 
  geom_text_repel(data=inv.labs,aes(label=name,y=-rr),size=4,alpha=1,point.padding=.5,min.segment.length = 0,nudge_y=-1.5,nudge_x=-.25) +
  facet_grid( ~ chrom,scales="free_x") + lims(y=c(-5,0)) +
  labs(x="",y="") + theme(panel.grid = element_blank(),axis.line=element_blank(),axis.ticks = element_blank(),axis.text=element_blank(),strip.text = element_blank(),
                          plot.margin=margin(t = 0, r = .5, b = .5, l = .5, unit = "pt"))

gg<-list(gg + theme(plot.margin=margin(t = .5, r = .5, b = 0, l = .5, unit = "pt")),gginv)
  }
  
  return(gg)
}


########################
draw_af_trajectories=function(afData,colorby="af0",rowcol="col",draw_cages=FALSE,mark_ts=TRUE,label_clNums=FALSE){
  ########################
  
  df <- afData %>% 
    filter(!is.na(tpt)) %>% mutate(startT=as.numeric(chop(as.character(timeseg)," ",2)),stopT=as.numeric(chop(as.character(timeseg)," ",4))) %>%
    group_by(chrom,pos,timeseg) %>%  
    mutate(af0=mean(af[tpt==0]),po=mean(af[tpt==startT]), pt=mean(af[tpt==stopT]), shift=pt - po) %>% 
    ungroup() %>% mutate(af=ifelse(shift<0,(af-af0)*-1,af-af0),
                         af0=ifelse(af0>0.5,1-af0,af0),
                         site=paste0(chrom,":",pos)
                         #s=calc_s_over_g(po,pt,3)
                         ) %>% 
    rename("col"=colorby)
  
  
  gg <- df %>%
    ggplot(aes(x=tpt,y=af))
  
  ## timeseg markers
  if(mark_ts){
    gg <- gg + geom_ribbon(data = df %>% filter(stopT-startT == 1,tpt==startT | tpt==stopT) %>% group_by(tpt,timeseg) %>% 
                             summarise(af=Inf) ,
                           aes(ymin=-Inf,ymax=Inf),alpha=.6,fill="azure3") 
  }
  
  ## per-cage lines
  if(draw_cages){
    gg<- gg +  geom_line(aes(group=paste0(site,cage)),alpha=.2,size=.5,color="gray") 
  }
  ## cross-cage average lines
  gg <- gg + geom_line(data=df %>% group_by(tpt,site,timeseg,col) %>% summarise(af=mean(af)),
                       aes(group=site,color=col),alpha=.2,size=1.5)
  
  ## highlight segment
  gg <- gg + geom_line(data = df %>% filter(stopT-startT == 1) %>% group_by(tpt,site,timeseg,startT,stopT,col) %>% 
                         summarise(af=mean(af))  %>% filter(tpt==startT | tpt==stopT),
                       aes(group=site,color=col),alpha=.6,size=1.5) 
  ## add cluster labels
  if(label_clNums){
    gg <- gg + geom_text_repel(data=df %>% filter(tpt==5) %>% 
                                 group_by(chrom,tpt,cl,timeseg) %>% 
                                 summarise(af=mean(af)) %>% 
                                 mutate(clNum=paste0(chop(cl,"\\.",1),".",chrom)),
                               aes(label=clNum),size=2 )
  }
  
  ## add color label
  if(colorby=="af0"){
    if(rowcol=="row"){ gg <-gg + labs(color="Initial Minor Allele Frequency")}
    if(rowcol=="col"){ gg <-gg + labs(color="Initial MAF") }
  }
  if(colorby=="s"){
    gg <-gg + 
      labs(color="Estimated\nSelection\nCoefficient") 
  }
  
  ## set orientation
  if(rowcol=="row"){gg <- gg + facet_wrap(~timeseg,nrow=1)}
  if(rowcol=="col"){gg <- gg + facet_wrap(~timeseg,ncol=1)}
  
  ## final fomatting
  gg <- gg +theme_minimal() + scale_color_viridis_c() +
    labs(x="Timepoint",y="Rising Allele Frequency (Difference from Baseline)")
  
  
  return(gg)
}
###################
