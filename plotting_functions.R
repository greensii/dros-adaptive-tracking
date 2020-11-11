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

plot_sites_manhattan=function(chrom,pos,score,group=NULL,ylab="score",mytitle="",pointSize=.7,pointShape=16,pointAlpha=.7){
  df=data.frame(chrom,pos,score);
  if(!is.null(group)){df$group=group}
  gg<-ggplot(df,aes(x=pos/1000000,y=score)) + facet_wrap(vars(chrom),scales = "free_x",nrow=1) + 
    theme_minimal() + labs(y=ylab,x="position (Mb)") + ggtitle(mytitle)
  if(!is.null(group)){
    gg<- gg + geom_point(aes(color=group),size=pointSize,shape=pointShape,alpha=pointAlpha)
  } else{
    gg<-gg + geom_point(size=pointSize,shape=pointShape,pointAlpha=.7)
  }
  gg
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