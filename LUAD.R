if (T) {
 
  dir.create("Files")
  dir.create("data")
  dir.create("Figures")
  dir.create("00_origin_datas/GEO",recursive = T)
  dir.create("00_origin_datas/TCGA")
  dir.create("00_pre_datas/GEO",recursive = T)
  dir.create("00_pre_datas/TCGA")
}
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)

get_riskscore<-function(dat,os,os.time,step=T,direction=c("both", "backward", "forward")[1]){
  tcga_dat1 <- cbind(time=os.time,
                     status=os,
                     dat)
  tcga_dat1=crbind2DataFrame(tcga_dat1)
  colnames(tcga_dat1)=gsub('-','__',colnames(tcga_dat1))
  gene111=gsub('-','__',colnames(dat))
  fmla <- as.formula(paste0("Surv(time, status) ~"
                            ,paste0(gene111,collapse = '+')))
  cox <- coxph(fmla, data =as.data.frame(tcga_dat1))
  if(step==T){
    cox1 <- step(cox,direction =direction)
  }else{
    cox1=cox
  }
  lan <- coef(cox1)
  #round(lan, 3)
  genes <- names(cox1$coefficients)
  mult_results=paste0(round(lan, 3), '*', names(lan),collapse = '+')
  risk.tcga <- as.numeric(lan%*%as.matrix(t(tcga_dat1[,genes])))
  
  data_gene_score_final<-tcga_dat1
  data_gene_score_final$Samples<-rownames(data_gene_score_final)
  data_gene_score_final$riskscore=risk.tcga
  data_gene_score_final$riskscorez=mosaic::zscore(risk.tcga)
  #
  optimalCutoff <- survminer::surv_cutpoint(data.frame(time=data_gene_score_final$time/365,
                                                       event=data_gene_score_final$status,
                                                       risk=data_gene_score_final$riskscore), 
                                            time = "time", event = "event",variables = c("risk"))
  optimalCutoff=optimalCutoff$cutpoint$cutpoint[1]
  #print(optimalCutoff)
  #optimalCutoff=median(data_gene_score_final$riskscore)
  #optimalCutoff=0
  data_gene_score_final$Risk=ifelse(data_gene_score_final$riskscore>optimalCutoff,'High','Low')
  table(data_gene_score_final$Risk)
  data_gene_score_final$cutoff=optimalCutoff
  return(list(result=data_gene_score_final,module.gene=genes,model=mult_results))
}
my_mutiviolin=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                       #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                       bw=T,xlab='',ylab='score',title='',size=3,angle = 45, hjust = 1,
                       legend.position='top',fill='group',notch=F){
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(dat.melt,aes(x=type, y=value,fill=Group)) +
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill='Group')+
    geom_violin(trim = F,position=position_dodge(0.8))+  
    scale_fill_manual(values = group_cols)+
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust))
  return(p)
}
plotMutiBar <-plotMutiBar <-function(dat=Age1_compare,ist=F,margin=T,xlb='',ylb='',lineCol='black',lineW=0.5,legTitle='Group',showValue=F,showLine=T,xangle=0,isAuto=T,color){
  library(ggplot2)
  #library(tidyverse)
  #library(reshape2)
  #library(optparse)
  if(ist){
    dat=t(dat)
  }
  lbc=colnames(dat)
  lbr=row.names(dat)
  bk_dat=dat
  if(margin){
    dat=dat%*%diag(1/c(apply(t(dat), 1, sum)))
  }
  row.names(dat)=paste0('R',1:(nrow(dat)))
  colnames(dat)=paste0('C',1:(ncol(dat)))
  row.names(bk_dat)=paste0('R',1:(nrow(bk_dat)))
  colnames(bk_dat)=paste0('C',1:(ncol(bk_dat)))
  #df=cbind(bg=paste0('R',1:nrow(dat)),dat)
  #colnames(df)=c('bg',paste0('C',1:(ncol(dat))))
  tp.dat=as.data.frame(cbind(bg=row.names(dat),dat))
  tp.dat[,1]=as.character(tp.dat[,1])
  for(i in 2:ncol(tp.dat)){
    tp.dat[,i]=as.numeric(as.character(tp.dat[,i]))
  }
  mt.df=reshape2::melt(tp.dat)
  colnames(mt.df)=c('bg','variable','value')
  
  pg=ggplot(mt.df, aes(x=variable, y=value, fill=bg))+geom_bar(stat = "identity", width=lineW, col=lineCol) 
  if(showLine){
    for (i in 2:(ncol(tp.dat)-1)) {
      tmp=tp.dat[order(tp.dat[,1],decreasing = T),]
      tmp[,i]=base::cumsum(tmp[,i])
      tmp[,i+1]=base::cumsum(tmp[,i+1])
      colnames(tmp)[c(i,i+1)]=c('STY','ED')
      tmp1=cbind(tmp,STX=rep(i-1+lineW/2,nrow(tmp))
                 ,EDX=rep(i-lineW/2,nrow(tmp)))
      pg=pg+geom_segment(data=tmp1,aes(x=STX, xend=EDX, y=STY, yend=ED))
    }
  }
  
  if(showValue){
    pg=pg+geom_text(data=mt.df,aes(label=sprintf("%0.2f", round(value, digits = 2))),position=position_stack(vjust=0.5))
  }
  pg=pg+scale_x_discrete(breaks = paste0('C',1:(ncol(dat))),label = lbc)
  pg=pg+labs(x=xlb, y=ylb)+ggsci::scale_fill_npg()+theme(legend.position = "bottom")
  pg=pg+ggsci::scale_fill_npg()+scale_fill_manual(values = color, 
                                                  breaks = paste0('R', 1:nrow(dat)), 
                                                  labels = lbr, 
                                                  name = legTitle) 
  if(xangle>0){
    pg=pg+theme(axis.text.x = element_text(angle = xangle, hjust = 1),legend.position = "bottom")
  }
  
  g.tb=matrix(0,nrow=ncol(dat),ncol=ncol(dat))
  for(i in 1:(ncol(dat))){
    for(j in 1:ncol(dat)){
      if(i!=j){
        g.tb[i,j]=round(-log10((chisq.test(bk_dat[,c(i,j)])$p.value)),2)
      }
    }
  }
  colnames(g.tb)=lbc
  row.names(g.tb)=lbc
  g.tb=reshape2::melt(g.tb) 
  colnames(g.tb)=c('A1','A2','A3')
  g.tb$A4=paste0(g.tb[,3],ifelse(g.tb[,3]>-log10(0.05),'(*)',''))
  stable.p=ggplot(g.tb, aes(A1, A2)) + geom_tile(aes(fill = A3),colour = "white") +xlab('')+ylab('')+ scale_fill_gradient(low = "white",high = "steelblue")+geom_text(aes(x=A1,y=A2,label=A4))+theme(legend.position="none",axis.title.x=element_blank(),axis.text.x=element_blank())
  stable.p=stable.p+ggtitle('-log10(CHl-Squared p value)')
  if(isAuto){
    g1=ggpubr::ggarrange(stable.p,pg, ncol = 1, nrow = 2,heights = c(0.5,1),align = "hv")
    return(g1)
  }else{
    return(list(Bar=pg,Table=stable.p))
  }
}

wb_beeswarm_plot <- function(dat = NULL,
                             show_compare = T,
                             xlab = 'Groups',
                             ylab = '',
                             method = c('t.test', 'wilcox.test')[1],
                             col = mycolor,
                             leg.pos = c('top','left','right','bottom','none')[1],
                             title = NULL,
                             group = 'Cluster') {
  library(ggbeeswarm)
  colnames(dat) <- c('Feature', 'Cluster')
  
  
  p1 <- ggplot(dat, aes(Cluster, Feature, color = Cluster)) + geom_quasirandom(method = "frowney") +
    ggtitle(title) + scale_color_manual(values = col[1:length(unique(dat$Cluster))]) +
    xlab(xlab) + ylab(ylab) + guides(color=guide_legend(title = group)) + theme_classic() +
    theme(legend.position=leg.pos)
  
  
  if(show_compare){
    uni.group = as.character(unique(dat$Cluster))
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,
                                     method = method,
                                     label= "p.signif", 
                                     step_increase = 0.0)
  }
  return(p1)
}
bioForest=function(rt=null,col){
  #
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  
  #
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  #
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  #
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin(trim = F)+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_boxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                    fill= "Group",label=c("p.format",'p.signif')[1],notch=F,
                    xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat=ARGs.score[tcga.subtype$Samples,'score']
  # group=tcga.subtype$Cluste
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data, aes(x=Group, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +
    scale_fill_manual(values = group_cols)+   #
    # if(length(names(table(group)))>2){
    #   test_method=''
    # }
    ggpubr::stat_compare_means(aes(group=Group), label = label, method = test_method)+
    labs(x=xlab, y = ylab, fill = fill) +
    theme_bw()+
    theme(legend.position =legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = x.size),
          axis.text.y = element_text(size = y.size)) # 
  return(p)
}

my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   #
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    # theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                 #
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) # 
  return(p)
}
my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","grey"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#
  return(p)
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  ## 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) ## 
  return(geneList)
}

coxFun <- function(dat){
  library(survival)
  colnames(dat)=c('time','status','gene')
  fmla=as.formula("Surv(time,status)~gene")
  cox=coxph(fmla,data=dat)
  p=summary(cox)[[7]][5]
  result=c(p,summary(cox)[[8]][1],summary(cox)[[8]][3],summary(cox)[[8]][4])
  return(result)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  set.seed(2021)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
wb_boxplot <- function(dat = data,
                       groups = groups,
                       xlab = '',
                       ylab = '',
                       xangle = 90,
                       title = 'Groups',
                       col = mycolor) {
  tmp.dat <- data.frame()
  for (ge in rownames(dat)) {
    print(ge)
    tmp <- data.frame(Samples = colnames(dat),
                      Genes = ge,
                      Values = as.numeric(dat[ge, ]),
                      Groups = groups)
    tmp.dat <- rbind(tmp.dat, tmp)
  }
  
  
  print(head(tmp.dat))
  library(ggpubr)
  if (length(unique(groups)) > 2) {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "anova", label = "p.signif") +
      scale_fill_manual(values = col) + theme_classic() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  } else {
    tmp_plot <- ggplot(tmp.dat, 
                       aes(x=Genes, y=Values, fill=Groups)) +
      geom_boxplot(notch = F) +  
      stat_compare_means(method = "t.test", label = "p.signif") +
      scale_fill_manual(values = col) + theme_classic() +
      theme(axis.text.x = element_text(angle=xangle, 
                                       hjust = 0.95,
                                       vjust = 0.95)) +
      xlab(xlab) + ylab(ylab) + guides(fill = guide_legend(title = title))
    return(tmp_plot)
  }
}

####TCGA-LUAD################
tcga.cli1<-read.delim('00_origin_datas/TCGA/Merge_LUAD_clinical.txt',sep='\t',header = T)
colnames(tcga.cli1)[1:20]
tcga.cli1$number_pack_years_smoked
table(tcga.cli1$tobacco_smoking_history)
tcga.cli1=data.frame(Samples=tcga.cli1$A0_Samples,
                     Age=tcga.cli1$A17_Age,
                     Gender=tcga.cli1$A18_Sex,
                     Smoking_history=tcga.cli1$tobacco_smoking_history,
                     T.stage=tcga.cli1$A3_T,
                     N.stage=tcga.cli1$A4_N,
                     M.stage=tcga.cli1$A5_M,
                     Stage=tcga.cli1$A6_Stage)
tcga.cli1$Samples=paste0(tcga.cli1$Samples,'-01')
rownames(tcga.cli1)=tcga.cli1$Samples
head(tcga.cli1)
table(tcga.cli1$Smoking_history)

table(tcga.cli1$T.stage)
tcga.cli1$T.stage=gsub('[ab]','',tcga.cli1$T.stage)
tcga.cli1$T.stage[tcga.cli1$T.stage=='TX']<-NA

table(tcga.cli1$N.stage)
tcga.cli1$N.stage[tcga.cli1$N.stage=='NX'|tcga.cli1$N.stage=='']<-NA

table(tcga.cli1$M.stage)
tcga.cli1$M.stage=gsub('[ab]','',tcga.cli1$M.stage)
tcga.cli1$M.stage[tcga.cli1$M.stage=='MX'|tcga.cli1$M.stage=='']<-NA

table(tcga.cli1$Stage)
tcga.cli1$Stage=gsub('[AB]','',tcga.cli1$Stage)
tcga.cli1$Stage[tcga.cli1$Stage=='']<-NA
tcga.cli1$Stage=gsub('Stage ','',tcga.cli1$Stage)


tcga.pancancer.cli=read.xlsx
head(tcga.pancancer.cli)
tcga.cli2=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='LUAD'),]
head(tcga.cli2)
tcga.cli2=data.frame(Samples=paste0(tcga.cli2$bcr_patient_barcode,'-01'),
                     tcga.cli2[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
head(tcga.cli2)
tcga.cli2$OS.time
tcga.cli2=tcga.cli2 %>% drop_na(OS.time)
####
tcga.cli2=tcga.cli2[tcga.cli2$OS.time>0,]
dim(tcga.cli2)

tcga.cli=merge(tcga.cli1,tcga.cli2,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
tcga.cli=as.data.frame(tcga.cli)
fivenum(as.numeric(tcga.cli$Age))
tcga.cli$Age1=ifelse(as.numeric(tcga.cli$Age)>66,'>66','<=66')
dim(tcga.cli)
# 509  17
head(tcga.cli)


#
tcga_data<-read.delim('00_origin_datas/TCGA/Merge_RNA_seq_FPKM _LUAD.txt',sep='\t',header = T,row.names = 1,check.names = F)
tcga_data[1:4,1:4]
table(substr(colnames(tcga_data),14,15))
tcga_data <- exp_ensg2symbol(tcga_data)

sample_T=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==1)]#
sample_N=colnames(tcga_data)[which(as.numeric(substr(colnames(tcga_data),14,15))==11)]#
length(sample_N)
length(sample_T)
tcga_type=data.frame(Samples=c(sample_T,sample_N),type=rep(c('Tumor','Normal'),c(length(sample_T),length(sample_N))))
rownames(tcga_type)=tcga_type$Samples
table(tcga_type$type)
# Normal  Tumor 
# 59    513

genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]


range(tcga_data)
tcga.exp.all=log2(tcga_data[intersect(rownames(tcga_data),mrna_genecode$SYMBOL),tcga_type$Samples]+1)
range(tcga.exp.all)
tcga.exp=tcga.exp.all[,intersect(tcga.cli$Samples,sample_T)]
dim(tcga.exp)
# 19503   500
tcga.cli=tcga.cli[intersect(tcga.cli$Samples,sample_T),]
dim(tcga.cli)
saveRDS(tcga.exp.all,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.all.RDS')
saveRDS(tcga.exp,file = '00_pre_datas/TCGA/LUAD_FPKM_tcga.exp.RDS')
########GSE31210#############
GSE31210 <- getGEOExpData('GSE31210')
saveRDS(GSE31210,file = "00_origin_datas/GEO/GSE31210.RDS")
#load('00_origin_datas/GEO/GSE31210.RData')

GSE31210.cli=GSE31210$Sample
GSE31210.cli=data.frame(Samples=GSE31210.cli$Acc,
                        Age=GSE31210.cli$`age (years)`,
                        Gender=GSE31210.cli$gender,
                        Status=GSE31210.cli$death,
                        OS.time=GSE31210.cli$`days before death/censor`)
rownames(GSE31210.cli)=GSE31210.cli$Samples
table(GSE31210.cli$Status)
GSE31210.cli=GSE31210.cli[which(GSE31210.cli$Status!='NULL'),]
GSE31210.cli$OS=ifelse(GSE31210.cli$Status=='alive',0,1)
GSE31210.cli <- GSE31210.cli[GSE31210.cli$OS.time>0,]
range(GSE31210.cli$OS.time)


GSE31210.exp=GSE31210$Exp$GPL570_54675_Data_col1
range(GSE31210.exp)
GSE31210.exp=log2(GSE31210.exp+1)
GSE31210.exp=exp_probe2symbol_v2(GSE31210.exp,GPL ='GPL570' )
range(GSE31210.exp)
dim(GSE31210.exp)
dim(GSE31210.cli)
GSE31210.exp=GSE31210.exp[,GSE31210.cli$Samples]
dim(GSE31210.exp)
#20549   226



####01.##########
dir.create('01_cluster')
ECS.gene <- read.table("01_cluster/EC.SENESCENCE.SIG.txt",sep = "\t")
ECS.gene <- ECS.gene$V1
length(ECS.gene)
#102
ECS.gene <- intersect(ECS.gene,rownames(tcga.exp))
######1.1#########
ECS.cox=cox_batch(tcga.exp[ECS.gene,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
ECS.cox[order(ECS.cox$p.value),]
table(ECS.cox$p.value<0.05)
#FALSE  TRUE 
#67    32
ECS.cox.fit=ECS.cox[ECS.cox$p.value<0.05,]
dim(ECS.cox.fit)
ECS.cox.fit$Type=ifelse(ECS.cox.fit$HR>1,'Risk','Portect')
sig_fit_gene=rownames(ECS.cox.fit)
length(sig_fit_gene)#32
writeMatrix(ECS.cox.fit,outpath = '01_cluster/sig_cox_res.txt')
writeMatrix(ECS.cox.fit,outpath = 'Files/sig_cox_res.txt')


########
pdf('01_cluster/sig_cox_bioForest.pdf',height = 12,width = 6)
bioForest(rt = ECS.cox.fit,col=c('#CC6600','#CC6699'))
dev.off()


#####1.2 #######################
library(ConsensusClusterPlus)
clusterAlg_name=c('hc','pam','km','kmdist')[1]
distance_name=c('pearson','spearman','euclidean','binary','maximum','canberra','minkowski')[2]
consen_gene=rownames(ECS.cox.fit)
length(consen_gene)#32
tcga_consen_data=as.matrix(tcga.exp[consen_gene,tcga.cli$Samples])
tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   #11,2  21,2
#tcga_consen_data=t(scale(t(tcga_consen_data),scale = F))   #
#tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, mean))#
# tcga_consen_data=sweep(tcga_consen_data,1,apply(tcga_consen_data, 1, median))  #
#tcga_consen_data=as.dist(1-cor(tcga_consen_data,method = 'pearson'))
tcga_consen_data=as.matrix(tcga_consen_data)
dim(tcga_consen_data)
tcga_clust_subtype <- ConsensusClusterPlus(tcga_consen_data
                                           , maxK = 10, reps = 500, pItem = 0.8
                                           , pFeature = 1
                                           , title = "TCGA_subtype"
                                           , clusterAlg = clusterAlg_name
                                           , distance = distance_name
                                           , plot = "pdf"
                                           , writeTable = T
                                           , seed = 123456)

k=2
cluster.color=pal_nejm()(8)[c(2,1)]

tcga.subtype <- data.frame(Samples = names(tcga_clust_subtype[[k]]$consensusClass),
                           Cluster=tcga_clust_subtype[[k]]$consensusClass)
tcga.subtype$Cluster=paste0('C',tcga.subtype$Cluster)
table(tcga.subtype$Cluster)
write.csv(tcga.subtype[order(tcga.subtype$Cluster),],'01_cluster/TCGA_subtype.csv',row.names = F)
tcga.subtype.cli=merge(tcga.subtype,tcga.cli,by='Samples')
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples


tcga.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                                        data = tcga.subtype.cli),
                           data=tcga.subtype.cli,
                           conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                           title='TCGA-LUAD',#ggtheme=custom_theme(),
                           linetype = c("solid", "dashed","strata")[1],
                           palette = cluster.color,
                           legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                           # legend = c(0.8,0.75), # 
                           legend.title = "")
tcga.subtype.km=mg_merge_plot(tcga.subtype.km$plot,tcga.subtype.km$table,ncol = 1,nrow = 2,heights = c(3,1),align = 'v')
ggsave("01_cluster/tcga.subtype.km.pdf",width = 6,height = 6)
tcga.subtype.km
#####带HR tcga.subtype.km2####

tcga.subtype.km2=ggplotKMCox(data.frame(tcga.subtype.cli$OS.time/365, 
                                        tcga.subtype.cli$OS,
                                        tcga.subtype.cli$Cluster),
                             
                             title='TCGA-LUAD',show_confint = T,
                             palette = cluster.color  )


"#0072B5FF" "#BC3C29FF"
ggsave("01_cluster/tcga.subtype.km2.pdf",tcga.subtype.km2,width = 6,height = 6)
######1.3GSE31210#######
GSE31210_conse<-GSE31210.exp[intersect(consen_gene,rownames(GSE31210.exp)),]
dim(GSE31210_conse)
GSE31210_conse=t(scale(t(GSE31210_conse),scale = F))

GSE31210_subtype <- ConsensusClusterPlus(GSE31210_conse
                                         , maxK = 10, reps = 500, pItem = 0.8
                                         , pFeature = 1
                                         , title = "GSE31320_subtype"
                                         , clusterAlg = clusterAlg_name
                                         , distance = distance_name
                                         # , innerLinkage = 'ward.D2'
                                         , plot = "pdf"
                                         , writeTable = T
                                         , seed = 123456)
GSE31210.subtype <- data.frame(Samples = names(GSE31210_subtype[[k]]$consensusClass),Cluster=GSE31210_subtype[[k]]$consensusClass)
GSE31210.subtype$Cluster=paste0('C',GSE31210.subtype$Cluster)
rownames(GSE31210.subtype)=GSE31210.subtype$Samples
GSE31210.subtype$Cluster=gsub('C1','IS2',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('C2','IS1',GSE31210.subtype$Cluster)
GSE31210.subtype$Cluster=gsub('IS','C',GSE31210.subtype$Cluster)
table(GSE31210.subtype$Cluster)
writeMatrix(GSE31210.subtype,row=T,header=T,outpath = '01_cluster/GSE31210.subtype.txt')
#GSE31210.subtype=readMatrix(inpath = 'GSE31210.subtype.txt',row=T,header=T)
#KM
GSE31210.subtype.cli <- cbind(GSE31210.cli, GSE31210.subtype[GSE31210.cli$Samples, ])
GSE31210.subtype.cli=crbind2DataFrame(GSE31210.subtype.cli)
write.table(GSE31210.subtype.cli,'01_cluster/GSE31210.subtype.cli.txt',quote = F,row.names = F,sep='\t')
GSE31210.subtype.km=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Cluster,
                                            data = GSE31210.subtype.cli),
                               data=GSE31210.subtype.cli,
                               conf.int = F,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                               title='GSE31210',#ggtheme=custom_theme(),
                               linetype = c("solid", "dashed","strata")[1],
                               palette = cluster.color,
                               legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                               # legend = c(0.8,0.75), # 
                               legend.title = "")
GSE31210.subtype.km=mg_merge_plot(GSE31210.subtype.km$plot,GSE31210.subtype.km$table,ncol = 1,nrow = 2,heights = c(3,1),align = 'v')
GSE31210.subtype.km
ggsave("01_cluster/GSE31210.subtype.km.pdf",width = 6,height = 6)
#########

GSE31210.subtype.km2=ggplotKMCox(data.frame(GSE31210.subtype.cli$OS.time/365, 
                                            GSE31210.subtype.cli$OS,
                                            GSE31210.subtype.cli$Cluster),
                                 
                                 title='GSE31210',show_confint = T,
                                 palette = cluster.color  )
ggsave("01_cluster/GSE31210.subtype.km2.pdf",GSE31210.subtype.km2,width = 6,height = 6)
#####TCGA ####
library(cluster)  # 
library(ggplot2)  # 

# 
k_values <- 2:10  # 
silhouette_scores <- list()

# 
for (k in k_values) {
  #k=3
  # 
  cluster_labels <- as.numeric(factor(tcga_clust_subtype[[k]]$consensusClass))
  
  # 
  #dist_matrix <- dist(tcga_consen_data, method = "euclidean")
  dist_matrix <- as.dist(1 - tcga_clust_subtype[[2]]$consensusMatrix)
  # 
  silhouette_res <- silhouette(cluster_labels, dist_matrix)
  
  # 
  silhouette_scores[[as.character(k)]] <- silhouette_res
}

#
silhouette_plot_data <- data.frame()

# 
for (k in k_values) {
  silhouette_res <- silhouette_scores[[as.character(k)]]
  temp_data <- data.frame(k = rep(k, nrow(silhouette_res)),
                          silhouette_score = silhouette_res[, 3])
  silhouette_plot_data <- rbind(silhouette_plot_data, temp_data)
}

# 
ggplot(silhouette_plot_data, aes(x = silhouette_score, fill = as.factor(k))) +
  geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
  labs(title = "Silhouette Score for Different K values",
       x = "Silhouette Score",
       y = "Frequency") +
  facet_wrap(~k, scales = "free") +
  theme_minimal()

# 
avg_silhouette_scores <- sapply(k_values, function(k) {
  silhouette_res <- silhouette_scores[[as.character(k)]]
  mean(silhouette_res[, 3])
})

# 
print(avg_silhouette_scores)
# 0.7085526 0.5336110 0.4866199 0.4079760 0.4049546 0.4108158 0.4199397 0.4048958 0.4344770

#########
library(cluster)  # 
library(ggplot2)  # 

# 

k_values <- 2:10  # 
silhouette_scores <- list()

# 
for (k in k_values) {
  #k=3
  # 
  cluster_labels <- as.numeric(factor(GSE31210_subtype[[k]]$consensusClass))
  
  # 
  #dist_matrix <- dist(tcga_consen_data, method = "euclidean")
  dist_matrix <- as.dist(1 - GSE31210_subtype[[2]]$consensusMatrix)
  # 
  silhouette_res <- silhouette(cluster_labels, dist_matrix)
  
  # 
  silhouette_scores[[as.character(k)]] <- silhouette_res
}

# 
silhouette_plot_data <- data.frame()

# 
for (k in k_values) {
  silhouette_res <- silhouette_scores[[as.character(k)]]
  temp_data <- data.frame(k = rep(k, nrow(silhouette_res)),
                          silhouette_score = silhouette_res[, 3])
  silhouette_plot_data <- rbind(silhouette_plot_data, temp_data)
}

# 
ggplot(silhouette_plot_data, aes(x = silhouette_score, fill = as.factor(k))) +
  geom_histogram(binwidth = 0.05, alpha = 0.6, position = "identity") +
  labs(title = "Silhouette Score for Different K values",
       x = "Silhouette Score",
       y = "Frequency") +
  facet_wrap(~k, scales = "free") +
  theme_minimal()

# 
avg_silhouette_scores <- sapply(k_values, function(k) {
  silhouette_res <- silhouette_scores[[as.character(k)]]
  mean(silhouette_res[, 3])
})

# 
print(avg_silhouette_scores)
#0.7884439 0.7213918 0.4394935 0.5383908 0.5076547 0.5274752 0.5436681 0.5326163 0.1661283

# 02 ######################
################
###
tcga.purity=readxl::read_xlsx('data/PMID_34019806.xlsx',sheet = 'Pan_TCGA')
tcga.purity=crbind2DataFrame(tcga.purity)
colnames(tcga.purity)[1]='Samples'
tcga.purity=tcga.purity[-1,]
tcga.purity$Samples=paste0(tcga.purity$Samples,'-01')
rownames(tcga.purity)=tcga.purity$Samples
head(tcga.purity)
table(tcga.purity$TCGA_project)
tcga.purity=tcga.purity[which(tcga.purity$TCGA_project=='LUAD'),]
dim(tcga.purity)

TMB.data <- data.frame(tcga.purity[tcga.subtype$Samples,'TMB'],tcga.subtype$Cluster)
pdf("02_cluster_SNV/fig2a.pdf",height = 6,width = 6)
wb_beeswarm_plot(dat = TMB.data,
                 show_compare = T,
                 xlab = 'Groups',
                 ylab = '',
                 method = c('t.test', 'wilcox.test')[1],
                 col = cluster.color,
                 leg.pos = c('top','left','right','bottom','none')[1],
                 title = NULL,
                 group = 'Cluster') 
dev.off()
#####
tcga.maf=getTCGAMAFByCode('LUAD')
tcga.subtype.use=tcga.subtype
table(tcga.subtype.use$Cluster)
colnames(tcga.subtype.use)[1]='Tumor_Sample_Barcode'
tcga.subtype.use$Tumor_Sample_Barcode=substr(tcga.subtype.use$Tumor_Sample_Barcode,1,12)
tcga.subtype.use.C1=tcga.subtype.use[which(tcga.subtype.use$Cluster=='C1'),]
tcga.subtype.use.C2=tcga.subtype.use[which(tcga.subtype.use$Cluster=='C2'),]
write.table(tcga.subtype.use.C1,file='02_cluster_SNV/tcga.subtype.c1.txt')
write.table(tcga.subtype.use.C2,file='02_cluster_SNV/tcga.subtype.c2.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.subtype.use.C1$Tumor_Sample_Barcode))
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = '02_cluster_SNV/tcga.subtype.c1.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.subtype.use.C2$Tumor_Sample_Barcode))
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = '02_cluster_SNV/tcga.subtype.c2.txt')
tcga.maf2@clinical.data



tcga.mut.dat <- tcga.maf
tcga.mut.dat <- as.data.frame(tcga.mut.dat@data)
tcga.mut.dat <- tcga.mut.dat[, c("Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification")]

tcga.mut.dat$Variant_Classification <- 1
tcga.mut.dat <- reshape2::dcast(data = tcga.mut.dat, Hugo_Symbol ~ Tumor_Sample_Barcode)
class(tcga.mut.dat)
rownames(tcga.mut.dat) <- tcga.mut.dat$Hugo_Symbol
tcga.mut.dat <- tcga.mut.dat[, -1]

colnames(tcga.mut.dat) <- paste0(colnames(tcga.mut.dat), '-01')
mut.samples <- intersect(colnames(tcga.mut.dat), tcga.subtype.cli$Samples)


tcga.mut.dat <- tcga.mut.dat[, mut.samples]
rownames(tcga.subtype.cli)=tcga.subtype.cli$Samples
tcga_mut_cli <- tcga.subtype.cli[mut.samples, ]

tcga.mut.dat.freq <- as.data.frame(rowSums(tcga.mut.dat))
colnames(tcga.mut.dat.freq) <- 'Freq'
tcga.mut.dat.freq$Genes <- rownames(tcga.mut.dat.freq)
library(dplyr)
str(tcga.mut.dat.freq)
head(tcga.mut.dat.freq)
tcga.mut.dat.freq <- dplyr::arrange(tcga.mut.dat.freq, desc(Freq))
head(tcga.mut.dat.freq)
dim(tcga.mut.dat.freq)

mut.genes <- rownames(tcga.mut.dat.freq)[tcga.mut.dat.freq$Freq > 3]
length(mut.genes)

tcga.mut.dat <- ifelse(tcga.mut.dat > 0, 'Mutant', 'WildType')
dim(tcga.mut.dat)

mut.res <- data.frame(C1 = NA,
                      C2 = NA
)
mut.p <- c()

for (ge in mut.genes) {
  #  print(ge)
  tmp <- table(tcga.mut.dat[ge, ], tcga_mut_cli$Cluster)
  pvalue <- fisher.test(tmp)
  mut.p <- c(mut.p, pvalue$p.value)
  mut.res <- rbind(mut.res, tmp[1, ])
}
mut.res <- na.omit(mut.res)
rownames(mut.res) <- mut.genes
class(mut.res)
mut.res$P.value <- mut.p

table(mut.res$P.value < 0.01)
# FALSE  TRUE 
# 9735   169
mut.res.filtered <- mut.res[which(mut.res$P.value < 0.01), ]
mut.res.filtered
dim(mut.res.filtered)
writeMatrix(mut.res.filtered,'02_cluster_SNV/subtype.mut.gene.txt')
cluster.color
pdf('02_cluster_SNV/Fig2b.pdf',height = 5,width =6)
oncoplot(maf=tcga.maf1,clinicalFeatures = 'Cluster',
         genes = rownames(mut.res.filtered)[1:15],
         sortByAnnotation = T,
         annotationColor = list(Cluster=c(C1=as.character(cluster.color[1]))))
dev.off()
pdf('02_cluster_SNV/Fig2c.pdf',height = 5,width =6)
oncoplot(maf=tcga.maf2,clinicalFeatures = 'Cluster',
         genes = rownames(mut.res.filtered)[1:15],
         sortByAnnotation = T,
         annotationColor = list(Cluster=c(C2=as.character(cluster.color[2]))))
dev.off()


pdf('02_cluster_SNV/Fig2d.pdf',height = 5,width =6)
OncogenicPathways(maf = tcga.maf1)
dev.off()

pdf('02_cluster_SNV/Fig2e.pdf',height = 5,width =6)
OncogenicPathways(maf = tcga.maf2)
dev.off()


####03_###########
dir.create('03_diff_genes')
tcga.limma=mg_limma_DEG(exp =tcga.exp[,tcga.subtype.cli$Samples],
                        group=tcga.subtype.cli$Cluster,ulab='C1',dlab = 'C2')
tcga.limma$Summary

#######
tcga.cluster.degs=tcga.limma$DEG
tcga.cluster.degs=tcga.cluster.degs[abs(tcga.cluster.degs$logFC)>log2(2) & tcga.cluster.degs$adj.P.Val<0.05,]
dim(tcga.cluster.degs)
write.csv(tcga.cluster.degs,'03_diff_genes/tcga.cluster.degs.csv')

######3.1####
p_cutoff <-0.05 
fc_cutoff <- log2(2)
degs_dat=tcga.limma$DEG
degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                            ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))

library(ggbreak)
library(ggplot2)
library(ggprism)

col=c("#FF6600","#3366CC","grey")
ylab='-log10 (adj.PVal)'
xlab='log2 (FoldChange)'
leg.pos='right'
plot3a<- ggplot(degs_dat, aes(x=logFC, y=-log10(adj.P.Val), color=type)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=col) +
  theme_bw() +
  theme(legend.position = leg.pos) +
  ylab(ylab) +
  xlab(xlab) +
  geom_vline(xintercept=c(-fc_cutoff,fc_cutoff), lty=3, col="black", lwd=0.5) +
  geom_hline(yintercept = -log10(p_cutoff), lty=3, col="black", lwd=0.5) 
  #coord_cartesian(ylim=c(0, 25)) # 
  #scale_y_break(c(25,100),#
                # space = 0.3,#
                # scales = 0.5)+#
  #theme_prism(palette = "black_and_white",
              # base_fontface = "plain", 
              # base_family = "serif", 
              # base_size = 16,
              # base_line_size = 0.8,
              # axis_text_angle = 0)
ggsave('03_diff_genes/fig3a.pdf',plot3a,height = 6,width = 6)

######3.2 ########
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.subtype.cli$Samples],group=tcga.subtype.cli$Cluster
                        ,ulab='C1',dlab = 'C2')

h.all.gmt<-read.gmt("data/h.all.v2023.1.Hs.entrez.gmt")
tcga.hallmark.gsea<-GSEA(tcga.geneList,TERM2GENE = h.all.gmt,seed=T)
library(enrichplot)
library(ggplot2)
gsea.p=clusterProfiler::dotplot(tcga.hallmark.gsea,split=".sign",showCategory=nrow(tcga.hallmark.gsea@result),
                                title='C1 vs C2',font.size=10)+facet_grid(~.sign)+
  theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
write.csv(tcga.hallmark.gsea@result,file = '03_diff_genes/cluster_GSEA.csv',quote = F)


ggsave('03_diff_genes/fig3b2.pdf',gsea.p,height = 8,width = 8)
#######04.########
dir.create('04_model')

sig.gene.cox=cox_batch(dat = tcga.exp[intersect(rownames(tcga.cluster.degs),rownames(tcga.exp)),tcga.cli$Samples],
                        time = tcga.cli$OS.time,event = tcga.cli$OS)
sig.gene.cox



table(sig.gene.cox$p.value<0.01)
# FALSE  TRUE 
# 98      105
pre.genes=rownames(sig.gene.cox[sig.gene.cox$p.value<0.01,])
length(pre.genes)#105
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
write.csv(sig.gene.cox,'04_model/sig.cox.csv')
######LASSO
tcga.lasso=get_riskscore.lasso(dat = tcga_model_data[,-c(1:2)],
                               os = tcga_model_data$OS,
                               os.time = tcga_model_data$OS.time)
length(tcga.lasso$lasso.gene)#11
tcga.lasso$plot

#####
fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(tcga.lasso$lasso.gene,collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)

lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
#"0.138*ANLN+-0.072*SLC34A2+0.097*ANGPTL4+0.113*FAM83A+0.118*GJB3"
write.csv(data.frame(gene=names(lan),coef=as.numeric(lan)),'04_model/gene_coef.csv',row.names = F)

####
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(lan), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))
survminer::ggforest(cox,data=tcga_model_data)
###########
lan.dataframe <- as.data.frame(lan)

lan.dataframe$gene <- rownames(lan.dataframe) 
lan.dataframe$gene <- factor(lan.dataframe$gene,levels = rownames(lan.dataframe)[order(lan.dataframe$lan)])
# 
lan.dataframe$color_group <- ifelse(lan.dataframe$lan > 0, "Positive", "Negative")
library(ggplot2)

# 


# 
p <- ggplot(lan.dataframe, aes(x=gene, y=lan,fill=color_group)) +
  geom_bar(stat="identity") +
  xlab("Gene Name") +
  ylab("Coefficient") +
  ggtitle("Gene Coefficients") +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#FF9999", "Negative" = "#6666CC")) +
  theme_bw()+
  guides(fill=FALSE)
p1 <- p+geom_text(aes(label=sprintf("%.3f", lan)), hjust=-0.2, size=3, color="black")
ggsave('04_model//gene_ Coefficients.pdf',height = 6,width = 9)


####
risktype.col=c('#CC99CC',"#669999")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
#######
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='TCGA-LUAD',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,
                      ylab='Overall Survival(OS)',
                      legend=c(0.85,0.8),#
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')

tcga.km.OS
######tcga.km.OS2 HR#####
tcga.km.OS2=ggplotKMCox(data.frame(tcga.risktype.cli$OS.time/365, 
                                   tcga.risktype.cli$OS,
                                   tcga.risktype.cli$Risktype),
                        
                        title='TCGA-LUAD',show_confint = T,
                        palette = risktype.col  )
###tcga#####
library(pROC)
library(ggplot2)

# 
roc_TCGA <- roc(tcga.risktype.cli$OS, tcga.risktype.cli$Riskscore)

# 
sensitivity_TCGA <- roc_TCGA$sensitivities
specificity_TCGA <- roc_TCGA$specificities

# 
best_threshold_TCGA <- roc_TCGA$thresholds[which.max(sensitivity_TCGA + specificity_TCGA)]

# 
predictions_TCGA <- ifelse(tcga.risktype.cli$Riskscore > best_threshold_TCGA, 1, 0)

# 
conf_matrix_TCGA <- table(Predicted = predictions_TCGA, Actual = tcga.risktype.cli$OS)

# 
TP_TCGA <- conf_matrix_TCGA[2, 2]  # 
FP_TCGA <- conf_matrix_TCGA[1, 2]  # 
TN_TCGA <- conf_matrix_TCGA[1, 1]  # 
FN_TCGA <- conf_matrix_TCGA[2, 1]  # 

# 
sensitivity_TCGA_value <- TP_TCGA / (TP_TCGA + FN_TCGA)
specificity_TCGA_value <- TN_TCGA / (TN_TCGA + FP_TCGA)
accuracy_TCGA_value <- (TP_TCGA + TN_TCGA) / sum(conf_matrix_TCGA)

# 
metrics_TCGA <- data.frame(
  Metric = c("Sensitivity", "Specificity", "Accuracy"),
  Value = c(sensitivity_TCGA_value, specificity_TCGA_value, accuracy_TCGA_value)
)

# 
print(metrics_TCGA)



#######
tcga.risktype.cli$Status=ifelse(tcga.risktype.cli$OS==0,'Alive','Dead')
tcga.model.p=my_riskplot(cli_dat = tcga.risktype.cli,cols =risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(tcga.risktype.cli$Riskscore),labs = '')


tcga.km.DSS=ggsurvplot(fit=survfit(Surv(DSS.time/365, DSS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Specific Survival(DSS)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DSS=mg_merge_plot(tcga.km.DSS$plot,tcga.km.DSS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DSS

tcga.km.DFI=ggsurvplot(fit=survfit(Surv(DFI.time/365, DFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',
                       title='TCGA-LUAD',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       ylab='Disease-Free interval(DFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.DFI=mg_merge_plot(tcga.km.DFI$plot,tcga.km.DFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DFI

tcga.km.PFI=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,fun = "pct",risk.table = T, 
                       size = 1,surv.median.line = 'hv',title='TCGA-LUAD',
                       legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       palette = risktype.col,
                       ylab='Progression-Free Interval(PFI)',
                       legend=c(0.85,0.8),#
                       ggtheme = theme_bw(base_size = 12) )
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.PFI


########
model.gene.df=data.frame(tcga.cli[,c('OS','OS.time')],t(tcga.exp[names(lan),tcga.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='TCGA',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.8,0.8),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=module.gene.km[[i]]$plot
}
tcga.module.km=mg_merge_plot(module.gene.km,ncol=3,nrow=2)



########4.2 ##########
GSE31210_model_data <- data.frame(GSE31210.cli[, c("OS.time", "OS")],
                             t(GSE31210.exp[intersect(names(lan),rownames(GSE31210.exp)), GSE31210.cli$Samples]))
colnames(GSE31210_model_data) <- gsub('-', '_', colnames(GSE31210_model_data))

risk.GSE31210=as.numeric(lan%*%as.matrix(t(GSE31210_model_data[GSE31210.cli$Samples,names(lan)])))

GSE31210.risktype.cli=data.frame(GSE31210.cli,Riskscore=risk.GSE31210)
GSE31210.risktype.cli$Risktype=ifelse(GSE31210.risktype.cli$Riskscore>median(risk.GSE31210),'High','Low')
GSE31210.roc=ggplotTimeROC(GSE31210.risktype.cli$OS.time,
                           GSE31210.risktype.cli$OS,
                           GSE31210.risktype.cli$Riskscore,mks = c(1,2,3))
GSE31210.roc
GSE31210.km=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = GSE31210.risktype.cli),
                       data=GSE31210.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE31210',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Overall Survival(OS)',
                       legend=c(0.85,0.25),#
                       ggtheme = theme_bw(base_size = 12))
GSE31210.km=mg_merge_plot(GSE31210.km$plot,GSE31210.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE31210.km
####HR GSE31210.km2######
GSE31210.km2=ggplotKMCox(data.frame(GSE31210.risktype.cli$OS.time/365, 
                                    GSE31210.risktype.cli$OS,
                                    GSE31210.risktype.cli$Risktype),
                         
                         title='GSE31210',show_confint = T,
                         palette = risktype.col  )
gse31210.model.p=my_riskplot(cli_dat = GSE31210.risktype.cli,cols = risktype.col,xlab = 'sample',
            a.ylab = 'Riskscore',b.labs = 'Time(days)',cutoff = median(GSE31210.risktype.cli$Riskscore),labs = '')



model.gene.df=data.frame(GSE31210.cli[,c('OS','OS.time')],t(GSE31210.exp[names(lan),GSE31210.cli$Samples]))
head(model.gene.df)
module.gene.km=list()
for (i in 1:length(names(lan))) {
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km[[i]]=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ group,
                                             data = model.gene.df1),
                                 data=model.gene.df1,
                                 conf.int = F,pval = T,risk.table = T, 
                                 fun = "pct",size = 1,surv.median.line = 'hv',
                                 title='GSE31210',legend.title=names(lan)[i],
                                 # legend.labs = c('High','Low'),
                                 linetype = c("solid", "dashed","strata")[1],
                                 palette = risktype.col,
                                 ylab='Overall Survival(OS)',
                                 legend=c(0.7,0.35),#
                                 ggtheme = theme_bw(base_size = 12))
  module.gene.km[[i]]=mg_merge_plot(module.gene.km[[i]]$plot,module.gene.km[[i]]$table,nrow=2,heights = c(3,1),align = 'v')
}
GSE31210.module.km=mg_merge_plot(module.gene.km,ncol=3,nrow=2)
ggsave("04_model/GSE31210.module.km.pdf",height =10 ,width = 10)
######HR GSE31210.module.km2####
module.gene.km2=list()
for (i in 1:length(names(lan))) {
  #i <- 1
  model.gene.df1=model.gene.df
  model.gene.df1$group=ifelse(model.gene.df1[,names(lan)[i]]>median(model.gene.df1[,names(lan)[i]]),'High','Low')
  module.gene.km2[[i]]=ggplotKMCox(data.frame(model.gene.df1$OS.time/365, 
                                              model.gene.df1$OS,
                                              model.gene.df1$group),
                                   
                                   title=names(lan)[i],show_confint = T,
                                   palette = risktype.col  )
}
GSE31210.module.km2=mg_merge_plot(module.gene.km2,ncol=3,nrow=2)
###GSE31210#####

library(pROC)
library(ggplot2)

# 
roc_GSE31210 <- roc(GSE31210.risktype.cli$OS, GSE31210.risktype.cli$Riskscore)

# 
sensitivity_GSE31210 <- roc_GSE31210$sensitivities
specificity_GSE31210 <- roc_GSE31210$specificities

# 
best_threshold_GSE31210 <- roc_GSE31210$thresholds[which.max(sensitivity_GSE31210 + specificity_GSE31210)]

# 
predictions_GSE31210 <- ifelse(GSE31210.risktype.cli$Riskscore > best_threshold_GSE31210, 1, 0)

# 
conf_matrix_GSE31210 <- table(Predicted = predictions_GSE31210, Actual = GSE31210.risktype.cli$OS)

# 
TP_GSE31210 <- conf_matrix_GSE31210[2, 2]  # 
FP_GSE31210 <- conf_matrix_GSE31210[1, 2]  # 
TN_GSE31210 <- conf_matrix_GSE31210[1, 1]  #
FN_GSE31210 <- conf_matrix_GSE31210[2, 1]  # 

# 
sensitivity_GSE31210_value <- TP_GSE31210 / (TP_GSE31210 + FN_GSE31210)
specificity_GSE31210_value <- TN_GSE31210 / (TN_GSE31210 + FP_GSE31210)
accuracy_GSE31210_value <- (TP_GSE31210 + TN_GSE31210) / sum(conf_matrix_GSE31210)

# 
metrics_GSE31210 <- data.frame(
  Metric = c("Sensitivity", "Specificity", "Accuracy"),
  Value = c(sensitivity_GSE31210_value, specificity_GSE31210_value, accuracy_GSE31210_value)
)

# 
print(metrics_GSE31210)

###########
fig4_1_HR=mg_merge_plot(mg_merge_plot(tcga.lasso$plot,p1,labels = c('','C')),
                   mg_merge_plot(tcga.model.p,tcga.roc,tcga.km.OS2,ncol=3,labels = LETTERS[4:6]),
                   #mg_merge_plot(tcga.km.DFI,tcga.km.DSS,tcga.km.PFI,ncol=3,labels = LETTERS[7:9]),
                   nrow=2)

fig4_2_HR=mg_merge_plot(mg_merge_plot(gse31210.model.p,GSE31210.roc,GSE31210.km2,ncol=3,labels = LETTERS[1:3]),
                     mg_merge_plot(module.gene.km2,ncol=3,nrow=2),nrow=2,heights = c(1,2),labels = c('','D'),align='hv')
ggsave('04_model/Fig4-1_HR.pdf',fig4_1_HR,height = 12,width = 20)
ggsave('04_model/figs4-2_HR.pdf',fig4_2_HR,height = 18,width = 20)


####05.##########
dir.create('05_Risktype.immune')
#  #############
dir.create('05_immune')
#######estimate####
tcga_estimate <- immu_estimate(exp = tcga.exp)

p5a <-  mg_PlotMutiBoxplot(tcga_estimate ,group = tcga.risktype.cli$Risktype,group_cols=risktype.col,legend.pos='top',test_method='wilcox.test',ylab = "score")



#####TIMER#####
tcga.timer=immu_timer(tcga.exp)
p5b<- mg_PlotMutiBoxplot(data = tcga.timer[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                         legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'score')

# mcp######
tcga.mcp <- immu_MCPcounter(exp = tcga.exp,isTCGA =T)
# mg_PlotMutiBoxplot(data = tcga.mcp[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
#                        legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'score')+labs(title = "MCP")


library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                    y=tcga.mcp[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',MCP_count=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
corrmat.color=colorRampPalette(c('blue', 'white','red'))(100)

p5c<-quickcor(tcga.mcp[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")
# ICGs####
tcga.icgs=immu_ICGs(tcga.exp)

icg.dat.RS=cbind(tcga.risktype.cli$Riskscore
                 ,tcga_model_data[tcga.risktype.cli$Samples,names(lan)]
                 ,tcga.icgs[tcga.risktype.cli$Samples,c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')])
#c('CTLA4','PDCD1','PDCD1LG2','LGALS9','CD80','CD28','HAVCR2')
colnames(icg.dat.RS)[1]='Riskcsore'

icg_cor_res <- Hmisc::rcorr(as.matrix(icg.dat.RS),type = 'spearman')
icg_cor_res$P[is.na(icg_cor_res$P)] <- 0
icg_cor_res.p=icg_cor_res$P
icg_cor_res.p[1:5,1:5]
icg_cor_res.p<-ifelse(icg_cor_res.p<0.0001,'****',
                      ifelse(icg_cor_res.p<0.001,'***', 
                             ifelse(icg_cor_res.p<0.01,'**',
                                    ifelse(icg_cor_res.p<0.05,'*',''))))

pdf('05_immune/p5d.pdf',height = 6,width = 7,onefile = F)
pheatmap(icg_cor_res$r[-c(1:6),c(names(lan),'Riskcsore')],
         color = circlize::colorRamp2(c(-1, 0, 1), c('#3B4992FF', 'white', '#EE0000FF')),
         main="Heatmap", # 
         display_numbers = icg_cor_res.p[-c(1:6),c(names(lan),'Riskcsore')], # 
         cluster_cols = F, # 
         cluster_rows = F,
         show_rownames = T, #
         show_colnames = T,
         fontsize_row = 12, # 
         fontsize_col = 16)
dev.off()

# TIDE  ###################
# tcga_tide_dat <- t(scale(t(tcga.exp),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = '05_immune/tcga_tide_dat.txt',quote = F, sep = '\t')
tcga_tide_res<-read.csv('05_immune/LIHC_0428.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
#tcga_tide_merge=cbind(tcga_tide_res[tcga.risktype.cli$Samples,],tcga.risktype.cli)
colnames(tcga_tide_dat)
colnames(tcga_tide_dat)
# tcga_tide_plot <-mg_PlotMutiBoxplot(tcga_tide_res[tcga.risktype.cli$Samples, c("TIDE",
#                                                                                        "Dysfunction",
#                                                                                        "Exclusion")],
#                                       ylab='Score',
#                               group = tcga.risktype.cli$Risktype,
#                                test_method='wilcox.test',
#                                legend.pos='top',
#                                group_cols = risktype.col)
# tcga_tide_plot
# ggsave(plot = tcga_tide_plot,
#        filename = '05_immune/tcga_tide_plot.pdf',
#        width = 6, height = 5)

#######
tcga_tide_res <- cbind(tcga_tide_res,
                       tcga.risktype.cli[rownames(tcga_tide_res), c("Riskscore", "Risktype")])


tcga_tide_cor <- cor_point(x = tcga_tide_res$TIDE,
                           y = tcga_tide_res$Riskscore,
                           xlab = 'TIDE',
                           ylab = 'RiskScore',top_col='#99CCFF',right_col='#99CC99'
)
tcga_tide_cor

tcga_Dysfunction_cor <- cor_point(x = tcga_tide_res$Dysfunction,
                                  y = tcga_tide_res$Riskscore,
                                  xlab = 'Dysfunction',
                                  ylab = 'RiskScore',
                                  top_col='#99CCFF',right_col='#99CC99')
tcga_Dysfunction_cor

tcga_Exclusion_cor <- cor_point(x = tcga_tide_res$Exclusion,
                                y = tcga_tide_res$Riskscore,
                                xlab = 'Exclusion',
                                ylab = 'RiskScore',
                                top_col='#99CCFF',right_col='#99CC99')
tcga_Exclusion_cor

# ######
p5ab<- cowplot::plot_grid(p5a,
                          p5b,
                          ncol = 2,labels = c("A",'B'),rel_widths = c(1,1.8),align ="h" )

ggsave('05_immune/p5ab.pdf',p5ab,height = 6,width = 15)

p5e <- cowplot::plot_grid(tcga_tide_cor,
                          tcga_Dysfunction_cor,
                          tcga_Exclusion_cor,
                          ncol = 3,labels = c("E",'F','G'))

ggsave('05_immune/p5e.pdf',p5e,height = 6,width = 15)

ggsave('05_immune/p5c.pdf',p5c,height = 6,width = 9)




####06.############
#######TMB##########
dir.create('06_mut')
tcga.tmb=mg_getTCGATMBByCode('LUAD')
tcga.tmb$Sample=paste0(tcga.tmb$Sample,'-01')
tcga.tmb.ri=merge(tcga.tmb,tcga.risktype.cli,by.x='Sample',by.y='Samples')
rownames(tcga.tmb.ri)=tcga.tmb.ri$Sample
head(tcga.tmb.ri)
fig6a=my_violin(dat = tcga.tmb.ri[-which(tcga.tmb.ri$TMB>10),'TMB'],
                group = tcga.tmb.ri$Risktype[-which(tcga.tmb.ri$TMB>10)],
                group_cols = risktype.col,
                fill ='Risktype',
                ylab = 'Tumor Mutation Burden',legend.position = 'none')
fig6a

#
tcga.tmb.cli=tcga.tmb.ri
library(survival)
library(survminer)
tmb.cutoff<-surv_cutpoint(thca.tmb.cli,
                          time="PFS.time",
                          event="PFS",
                          variables=c("TMB"))
summary(tmb.cutoff)
thca.tmb.cli$type <- ifelse(thca.tmb.cli$TMB > tmb.cutoff$cutpoint$cutpoint, 'TMB-High', 'TMB-Low')
tcga.tmb.cli$ri_tmb=rep('none',nrow(tcga.tmb.cli))
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-High')]='H-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-Low')]='L-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='Low' & tcga.tmb.cli$type=='TMB-High')]='L-Risk & H-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$Risktype=='High' & tcga.tmb.cli$type=='TMB-Low')]='H-Risk & L-TMB'
tcga.tmb.cli$ri_tmb[which(tcga.tmb.cli$ri_tmb=='none')]=NA
table(tcga.tmb.cli$ri_tmb)
fig6b=ggsurvplot(fit=survfit(Surv(OS.time, OS) ~ type,
                             data = data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                               , OS = tcga.tmb.cli$OS
                                               , type=tcga.tmb.cli$ri_tmb)),
                 data=data.frame(OS.time = tcga.tmb.cli$OS.time/365
                                 , OS = tcga.tmb.cli$OS
                                 , type=tcga.tmb.cli$ri_tmb),
                 conf.int = T,pval = T,fun = "pct",risk.table = T, size = 0.7,surv.median.line = 'hv',
                 title='TMB & Risktype',ggtheme=theme_classic(),
                 linetype = c("solid", "dashed","strata")[1],
                 palette = pal_lancet()(9)[2:5],
                 #legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                 legend = c(0.8,0.75), # 
                 legend.title = "")#,legend.labs =  c('TMB-High','TMB-Low')
fig6b
tmb.km.OS=mg_merge_plot(fig6b$plot,fig6b$table,nrow=2,heights = c(3,1),align = 'v')


##########
fig6ab=mg_merge_plot(fig6a,tmb.km.OS,
                     ncol=2,labels = c('A','B'),widths = c(1,1.3))


savePDF('06_mut/Fig6ab.pdf',fig6ab,height = 6,width = 12)
################
library(progeny)
tcga.pathway.activ=progeny(as.matrix(tcga.exp),scale = T)
dim(tcga.pathway.activ)
range(tcga.pathway.activ)
mg_PlotMutiBoxplot(tcga.pathway.activ[tcga.risktype.cli$Samples,]
                   , group = tcga.risktype.cli$Risktype
                   , legend.pos = 'top'
                   #, group_cols = cluster.color
                   , test_method =  c('kruskal.test','wilcox.test','anova')[2]
                   , add = 'boxplot'
                   , ylab = 'score')


pathway_cor_RS=cbind.data.frame(Riskscore=tcga.risktype.cli$Riskscore,
                                tcga.pathway.activ[tcga.risktype.cli$Samples,])
cor_res <- Hmisc::rcorr(as.matrix(pathway_cor_RS),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0

pdf('06_mut/pathway_cor_RS.pdf',height = 7,width = 7)
fig6c <- corrplot(as.matrix(cor_res$r),
                  p.mat = as.matrix(cor_res$P),
                  mar = c(0,0,1,1),diag = F,
                  col=COL2('PRGn'),#diverging_hcl(100,palette = 'Green-Orange')
                  tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
                  cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
                  cl.ratio = 0.1,cl.cex = 0.8,
                  addgrid.col = 'white',
                  method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
                  insig = 'label_sig',
                  sig.level=c(0.001,0.01,0.05),
                  pch.cex=1,is.corr=T,xpd=T)
dev.off()


##########################
tcga.risktype.cli.use=tcga.risktype.cli[,c('Samples','Risktype')]
colnames(tcga.risktype.cli.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.cli.use$Tumor_Sample_Barcode=substr(tcga.risktype.cli.use$Tumor_Sample_Barcode,1,12)
rownames(tcga.risktype.cli.use)=tcga.risktype.cli.use$Tumor_Sample_Barcode
tcga.risktype.cli.use=tcga.risktype.cli.use[order(tcga.risktype.cli.use$Risktype),]
write.table(tcga.risktype.cli.use,file='06_mut/tcga.risktype.cli.use.txt')

tcga.risktype.cli.use.h=tcga.risktype.cli.use[which(tcga.risktype.cli.use$Risktype=='High'),]
tcga.risktype.cli.use.l=tcga.risktype.cli.use[which(tcga.risktype.cli.use$Risktype=='Low'),]
write.table(tcga.risktype.cli.use.h,file='06_mut/tcga.Risktype.h.txt')
write.table(tcga.risktype.cli.use.l,file='06_mut/tcga.Risktype.l.txt')

tcga.maf=getTCGAMAFByCode('LUAD')
tcga.tmb.high=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,
                                               tcga.risktype.cli.use.h$Tumor_Sample_Barcode))
tcga.tmb.high<-read.maf(tcga.tmb.high@data,isTCGA=T,clinicalData = '06_mut/tcga.Risktype.h.txt')
tcga.tmb.high@clinical.data

tcga.tmb.low=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,
                                              tcga.risktype.cli.use.l$Tumor_Sample_Barcode))
tcga.tmb.low<-read.maf(tcga.tmb.low@data,isTCGA=T,clinicalData = '06_mut/tcga.Risktype.l.txt')
tcga.tmb.low@clinical.data

tcga.tmb_1=crbind2DataFrame(tcga.maf@data)
tcga.tmb.data=tcga.tmb_1[which(tcga.tmb_1$Tumor_Sample_Barcode %in% rownames(tcga.risktype.cli.use)),]
length(table(tcga.tmb.data$Tumor_Sample_Barcode))
#
tcga.tmb.data.uni.gene=unique(cbind(tcga.tmb.data$Hugo_Symbol,tcga.tmb.data$Tumor_Sample_Barcode))
tcga.tmb.data.uni.gene<-crbind2DataFrame(tcga.tmb.data.uni.gene)
colnames(tcga.tmb.data.uni.gene)<-c('GeneSymbol','Samples')
head(tcga.tmb.data.uni.gene)
#
gene.samp.cnt.tab=table(tcga.tmb.data.uni.gene[,2])
head(gene.samp.cnt.tab)


#
subtype_mut_num=data.frame(subtype=tcga.risktype.cli.use[intersect(names(gene.samp.cnt.tab),rownames(tcga.risktype.cli.use)),2]
                           ,`Number of Mut Genes` = as.numeric(gene.samp.cnt.tab))
head(subtype_mut_num)

#
colnames(tcga.risktype.cli.use)[1]='Samples'
tcga.tmb.data.uni.gene.1<-merge(tcga.tmb.data.uni.gene,tcga.risktype.cli.use,by.x = 'Samples',by.y = 'Samples')
gene.subty.cnt<-dcast(tcga.tmb.data.uni.gene.1,GeneSymbol ~ Risktype,length)
rownames(gene.subty.cnt)=gene.subty.cnt$GeneSymbol
gene.subty.cnt=gene.subty.cnt[,-1]
head(gene.subty.cnt)
#
gene.subty.cnt.s1=gene.subty.cnt[apply(gene.subty.cnt, 1,sum)>3,]
dim(gene.subty.cnt.s1)
#9780
head(gene.subty.cnt.s1)
##
s.mut.num=length(unique(tcga.tmb.data.uni.gene[,2]))
s.mut.num#495
#
subtype.mut.count=as.numeric(table(tcga.risktype.cli.use[intersect(unique(tcga.tmb.data.uni.gene[,2]),row.names(tcga.risktype.cli.use)),2]))
subtype.mut.count
# 249 246
pval.s.all=apply(gene.subty.cnt.s1, 1, function(x){
  pval.s=c()
  for(i in 1:dim(gene.subty.cnt.s1)[2]){
    pval=fisher.test(matrix(c(x[i]
                              ,sum(x)-x[i]
                              ,subtype.mut.count[i]-x[i]
                              ,s.mut.num-subtype.mut.count[i]-(sum(x)-x[i]))
                            ,nrow = 2,ncol = 2)
                     ,alternative = "greater")$p.value
    
    pval.s=c(pval.s,pval)
  }
  return(pval.s)
})
pval.s.all=t(pval.s.all)
head(pval.s.all)
#
gene.subty.cnt.s1.pval=cbind(gene.subty.cnt.s1,pval.s.all)
colnames(gene.subty.cnt.s1.pval)=c(colnames(gene.subty.cnt.s1),paste0(colnames(gene.subty.cnt.s1),'.pvalue'))
head(gene.subty.cnt.s1.pval)
gene.subty.cnt.s1.pval.sig=gene.subty.cnt.s1.pval[which(apply(pval.s.all, 1,function(x){return(sum(x<0.05))})>0),]
length(rownames(gene.subty.cnt.s1.pval.sig))#1091
top_gene=rownames(gene.subty.cnt.s1.pval.sig)[1:20]

top_gene_order <- as.data.frame(table(tcga.tmb.data.uni.gene[,1]))
top_gene_order <- top_gene_order[top_gene_order$Var1 %in%  top_gene,]
top_gene_order <- top_gene_order[order(top_gene_order$Freq,decreasing = T),]
top_gene <- top_gene_order$Var1
pdf('06_mut/Fig6d.pdf',height = 5,width = 10)
coOncoplot(m1=tcga.tmb.high, 
           m2=tcga.tmb.low, 
           m1Name="High",
           m2Name="Low",
           genes =  top_gene)
dev.off()


####07############
dir.create('07_immune_treatment')

#####IMvigor210#####
library("IMvigor210CoreBiologies")
data(cds)
pheno<-pData(cds)
head(pheno)
exper_tpm=mg_get_immu_pd1_treament_exp()
exper_id=exper_tpm$fpkm
exper_id$symbol=rownames(exper_id)
rownames(exper_id)<-exper_id$symbol
exper_id$symbol<-NULL
range(exper_id)
exper_id_use<-log2(exper_id+1)
dim(exper_id_use)
# 31085   348
range(exper_id_use)
#rownames(exper_id_use)=gsub('-','__',rownames(exper_id_use))
exper_id_use[1:5,1:5]
exper_id_use=exper_id_use[,rownames(pheno[which(pheno$binaryResponse!='NA'),])]
dim(exper_id_use)
# 31085   298
pheno=pheno[which(pheno$binaryResponse!='NA'),]

IMvigor210_model_data=data.frame(OS=pheno$censOS,OS.time = pheno$os,
                                 t(exper_id_use[intersect(names(lan),rownames(exper_id_use)),rownames(pheno)]))
head(IMvigor210_model_data)

fmla.IMvigor210 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                     ,paste0(names(lan),collapse = '+')))
cox.IMvigor210 <- coxph(fmla.IMvigor210, data =as.data.frame(IMvigor210_model_data))
IMvigor210_lan <- coef(cox.IMvigor210)

#IMvigor210_lan <- lan
risk.imv210=as.numeric(IMvigor210_lan%*%as.matrix(t(IMvigor210_model_data[rownames(pheno),names(IMvigor210_lan)])))

imv210.risktype.cli=cbind.data.frame(pheno,Riskscore=risk.imv210)
#######
imv210.data.point <- surv_cutpoint(imv210.risktype.cli, time = "os", event = "censOS",
                                   variables = 'Riskscore')
imv210.cutoff <- as.numeric(summary(imv210.data.point)[1])
imv210.cutoff
imv210.risktype.cli$Risktype=ifelse(risk.imv210>imv210.cutoff,'High','Low')

imv210.roc=ggplotTimeROC(imv210.risktype.cli$os,
                         imv210.risktype.cli$censOS,
                         imv210.risktype.cli$Riskscore,mks = c(1,2,3,4,5))
imv210.roc


imv210.km=ggsurvplot(fit=survfit(Surv(os, censOS) ~ Risktype,
                                 data = imv210.risktype.cli),
                     data=imv210.risktype.cli,
                     conf.int = F,pval = T,risk.table = T, 
                     fun = "pct",size = 1,surv.median.line = 'hv',
                     title='IMvigor210',legend.title='Risktype',
                     legend.labs = c('High','Low'),
                     linetype = c("solid", "dashed","strata")[1],
                     palette = risktype.col,ylab='Overall Survival(OS)',
                     legend=c(0.85,0.8),#
                     ggtheme = theme_classic())
p7a=mg_merge_plot(imv210.km$plot,imv210.km$table,nrow=2,heights = c(3,1),align = 'v')

#####HR p7a_HR#######
imv210.km2=ggplotKMCox(data.frame(imv210.risktype.cli$os, 
                                  imv210.risktype.cli$censOS,
                                  imv210.risktype.cli$Risktype),
                       
                       title='IMvigor210',show_confint = T,
                       palette = risktype.col  )




p7a_HR=imv210.km2

# head(imv210.risktype.cli)
# table(imv210.risktype.cli$binaryResponse)
# imv210.bar=plotMutiBar(table(imv210.risktype.cli$binaryResponse,imv210.risktype.cli$Risktype))
# imv210.bar

p7b=ggplot(imv210.risktype.cli,aes(x=`Best Confirmed Overall Response`,y=Riskscore,fill=`Best Confirmed Overall Response`))+
  geom_boxplot()+stat_compare_means(aes(group=`Best Confirmed Overall Response`), label = 'p.format', method = 'kruskal.test')+
  theme(legend.position = 'none',text = element_text(family = 'Times',size = 12))


#####GSE78220###############
library(stringr)
GSE78220_cli <- getGEOSampleData('GSE78220')
GSE78220_cli1 <- GSE78220_cli[, c("Acc", "Title", "anti-pd-1 response", 
                                  "overall survival (days)", "vital status")]
colnames(GSE78220_cli1) <- c(c("Samples", "Title", "Rresponse", 
                               "OS.time", "OS"))
GSE78220_cli1 <- na.omit(GSE78220_cli1)
table(GSE78220_cli1$OS)
GSE78220_cli1$OS <- ifelse(GSE78220_cli1$OS == 'Alive', 0, 1)
rownames(GSE78220_cli1) <- GSE78220_cli1$Title

GSE78220_exp <- openxlsx::read.xlsx('00_origin_datas/GEO/GSE78220/GSE78220_PatientFPKM.xlsx',
                                    sheet = 1)
rownames(GSE78220_exp) <- GSE78220_exp$Gene
GSE78220_exp <- GSE78220_exp[, -1]
colnames(GSE78220_exp) <- str_split_fixed(colnames(GSE78220_exp), '\\.', 3)[, 1]
boxplot(GSE78220_exp[, 1:5])

GSE78220_exp <- log2(GSE78220_exp + 1)
boxplot(GSE78220_exp[, 1:5])

GSE78220_model_data <- cbind(GSE78220_cli1[, c("OS.time", 'OS')],
                             t(GSE78220_exp[, rownames(GSE78220_cli1)]))
GSE78220.genes <- intersect(names(lan), colnames(GSE78220_model_data))
GSE78220.genes

GSE78220_model_data=GSE78220_model_data[,c('OS.time','OS',GSE78220.genes)]

# GSE78220.genes=gsub('-','__',GSE78220.genes)
# colnames(GSE78220_model_data)=gsub('-','__',colnames(GSE78220_model_data))

fmla.GSE78220 <- as.formula(paste0("Surv(OS.time, OS) ~"
                                   ,paste0(GSE78220.genes,collapse = '+')))
cox.GSE78220 <- coxph(fmla.GSE78220, data =as.data.frame(GSE78220_model_data))
GSE78220_lan <- coef(cox.GSE78220)

risk.GSE78220=as.numeric(GSE78220_lan%*%as.matrix(t(GSE78220_model_data[GSE78220_cli1$Title,names(GSE78220_lan)])))

GSE78220_model_data$RS <- risk.GSE78220
GSE78220.data.point <- surv_cutpoint(GSE78220_model_data, time = "OS.time", event = "OS",
                                     variables = 'RS')
GSE78220.cutoff <- as.numeric(summary(GSE78220.data.point)[1])
GSE78220.cutoff

GSE78220.roc <- ggplotTimeROC(GSE78220_model_data$OS.time / 365,
                              GSE78220_model_data$OS,
                              risk.GSE78220,
                              mks = c(1,2,2.5))


GSE78220.risktype.cli=data.frame(GSE78220_cli1,
                                 Riskscore=risk.GSE78220,
                                 Risktype=ifelse(risk.GSE78220>=GSE78220.cutoff,'High','Low'))
GSE78220.km <-ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                     data = GSE78220.risktype.cli),
                         data=GSE78220.risktype.cli,
                         conf.int = F,pval = T,risk.table = T, 
                         fun = "pct",size = 1,surv.median.line = 'hv',
                         title='GSE78220',legend.title='Risktype',
                         legend.labs = c('High','Low'),
                         linetype = c("solid", "dashed","strata")[1],
                         palette = risktype.col,ylab='Overall Survival(OS)',
                         legend=c(0.85,0.85),#
                         ggtheme = theme_classic()) 

p7c=mg_merge_plot(GSE78220.km$plot,GSE78220.km$table,nrow=2,heights = c(3,1),align = 'v')

####HR p7c_HR#####
GSE78220.km2=ggplotKMCox(data.frame(GSE78220.risktype.cli$OS.time/365, 
                                    GSE78220.risktype.cli$OS,
                                    GSE78220.risktype.cli$Risktype),
                         
                         title='GSE78220',show_confint = T,
                         palette = risktype.col  )




p7c_HR=GSE78220.km2



p7d <-ggplot(GSE78220.risktype.cli,aes(x=Rresponse,y=Riskscore,fill=Rresponse))+
  geom_violin(trim = F)+  
  scale_fill_manual(values = pal_nejm()(10)[3:5] )+
  geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
  stat_compare_means(aes(group=Rresponse), label = 'p.format', method = 'kruskal.test')+
  theme(legend.position = 'none',text = element_text(family = 'Times',size = 12))+xlab('Response')+
  theme_classic()



#############
library(pRRophetic)
library(ggplot2)
############### Cisplatin,
set.seed(12345)
predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga.exp)
                                              , "Cisplatin"
                                              , selection=1
                                              ,dataset = "cgp2016")
predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)

tcga_durg_ic50_res <- predictedPtype_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
length(drugs)
for (drug in drugs) {
  print(drug)
  set.seed(12345)
  tmpic50 <- pRRopheticPredict(as.matrix(tcga.exp)
                               , drug
                               , selection=1
                               , dataset = "cgp2016")
  tmpic50 <- data.frame(tmpic50)
  colnames(tmpic50) <- drug
  tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
}
# save(tcga_durg_ic50_res,file='07_immune_treatment/tcga_durg_ic50.RData')
# readRDS('07_immune_treatment/tcga_durg_ic50.RData')
head(tcga_durg_ic50_res)
tcga_durg_ic50_res=tcga_durg_ic50_res[,-1]

tcga.ic50.dat=cbind(tcga.risktype.cli[,'Riskscore'],tcga_durg_ic50_res[tcga.risktype.cli$Samples,])
head(tcga.ic50.dat)
colnames(tcga.ic50.dat)[1]='Riskscore'
ic50.cor.RS=Hmisc::rcorr(as.matrix(tcga.ic50.dat),type = 'spearman')
ic50.cor.RS.res=data.frame(Names=names(ic50.cor.RS$r['Riskscore',]),
                           cor=as.numeric(ic50.cor.RS$r['Riskscore',]),
                           p.val=as.numeric(ic50.cor.RS$P['Riskscore',]))
ic50.cor.RS.res=ic50.cor.RS.res[-1,]
head(ic50.cor.RS.res)
colnames(ic50.cor.RS.res)=c('Drugs','cor','pvalue')
ic50.cor.RS.res$Drugs=factor(ic50.cor.RS.res$Drugs,
                             levels = ic50.cor.RS.res$Drugs[order(ic50.cor.RS.res$cor,decreasing = T)], ordered=TRUE)
ic50.cor.RS.res$pvalue=ifelse(ic50.cor.RS.res$pvalue==0,1e-16,ic50.cor.RS.res$pvalue)
head(ic50.cor.RS.res)
rownames(ic50.cor.RS.res)=ic50.cor.RS.res$Drugs


fit.drug=as.character(ic50.cor.RS.res$Drugs[which(ic50.cor.RS.res$pvalue<0.05 & abs(ic50.cor.RS.res$cor) >0.3)])
length(fit.drug)

ic50.cor.RS.res$group=rep('no sign',nrow(ic50.cor.RS.res))
ic50.cor.RS.res$group[which(ic50.cor.RS.res$cor<0)]='Resistance'
ic50.cor.RS.res$group[which(ic50.cor.RS.res$cor>0)]='Sentivity'

p7e=ggplot(ic50.cor.RS.res[fit.drug,],aes(x=cor,y=Drugs,fill=group))+
  # geom_bar()+
  scale_fill_manual(values = c('#66CCCC','#FF99CC'))+
  xlab('Correlation coefficient with IC50')+ylab('Drugs')+
  geom_bar(stat="identity",position = "dodge") +
  theme(axis.text.y = element_text(size =10),
        axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = 'none')+coord_flip()
p7e
#####
p7e = ggplot(ic50.cor.RS.res[fit.drug,], aes(x = cor, y = Drugs, fill = pvalue)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_gradient(low = "#78B3CE", high = "#FF2929", name = "p-value") + # 
  xlab("Correlation coefficient with IC50") + 
  ylab("Drugs") +
  geom_text(aes(label = round(cor, 2)), hjust = 0.5, size = 5) + # 
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'right') +
  coord_flip() +
  theme_minimal() # 
p7e

fig7_HR=mg_merge_plot(mg_merge_plot(p7a_HR,p7b,ncol=2,widths = c(1,1.5),labels = LETTERS[1:2]),
                   mg_merge_plot(p7c_HR,p7d,ncol=2,widths = c(1,1.5),labels = LETTERS[3:4]),
                   p7e,labels = c('','','E'),heights = c(1,0.9,1),
                   nrow=3,ncol=1)
fig7_HR
savePDF('07_immune_treatment/Fig7_HR.pdf',fig7_HR,height = 15,width = 12)



save.image("LUAD4.Rdata")

###########
lan <- lan[-3]



####
risktype.col=c('#CC99CC',"#669999")

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.subtype.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.subtype.cli,Riskscore=risk.tcga)
#######
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc=ggplotTimeROC(tcga.risktype.cli$OS.time,
                       tcga.risktype.cli$OS,
                       tcga.risktype.cli$Riskscore,mks = c(1,3,5))
tcga.roc
tcga.km.OS2=ggplotKMCox(data.frame(tcga.risktype.cli$OS.time/365, 
                                   tcga.risktype.cli$OS,
                                   tcga.risktype.cli$Risktype),
                        
                        title='TCGA-LUAD',show_confint = T,
                        palette = risktype.col  )

tcga.km.OS2

