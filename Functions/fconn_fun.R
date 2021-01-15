require(ggplot2)
require(fields)
require(gplots)
require(reshape2)
require(dplyr)
require(circlize)

net2vec = function(innet) innet[lower.tri(innet)] 

subgatheraggfun = function(datain=subj.cor.df,indices1,indices2,positive=TRUE) {
  # if positive=TRUE, then average the correlations across all participants and all AF
  # and subset to positive edges. 
  
  subsetindices = datain$node2%in%indices1 & datain$node1%in%indices2
  
  submeancor = datain[subsetindices,]  
  
  # grab the other ones: as a result, node1=list1,node2=list2 produces same results as node1=list2,node2=list1
  subsetindices2 = datain$node2%in%indices2 & datain$node1%in%indices1
  submeancor = rbind(submeancor,datain[subsetindices2,])
  submeancor = submeancor[!duplicated(submeancor),]
  
  submeancorpos = submeancor
  
  if(positive) {
    submeancorpos$mean_acrossmb = rowMeans(submeancorpos[,c("SB.3.3.mm","SB.2.mm","MB.2","MB.3","MB.4","MB.6","MB.8","MB.9","MB.12")],na.rm=TRUE)
    submeancorpos$node12 = paste(submeancorpos$node1,submeancorpos$node2)
    submeancorpos = submeancorpos %>% group_by(node12) %>% mutate(mean_mean_cor = mean(mean_acrossmb,na.rm=TRUE))
    tempindices = submeancorpos$mean_mean_cor<0
    tempindices[is.na(tempindices)]=FALSE
    for (j in c("SB.3.3.mm","SB.2.mm","MB.2","MB.3","MB.4","MB.6","MB.8","MB.9","MB.12")) {}
    submeancorpos[tempindices,j] = NA
  }
  
  #tidyr retired gather; now suggests "pivot_longer"
  #subgather = gather(submeancorpos,key='MB',value='Corr',varnames)
  
  subgather = submeancorpos %>% pivot_longer(cols=c("SB.3.3.mm","SB.2.mm","MB.2","MB.3","MB.4","MB.6","MB.8","MB.9","MB.12"), names_to = "MB", values_to = "Corr")
  subgather$MB = factor(subgather$MB)
  subgather$MB = relevel(subgather$MB, ref="SB.2.mm")
  #subgather$Zcorr = tanh(subgather$Corr)
  
  # average the edges: 
  subgatheragg = aggregate(Corr~MB+id+gender,data=subgather,FUN=mean,na.rm=TRUE)
  subgatheragg$MB = factor(subgatheragg$MB,levels=c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12"))
  subgatheragg
}


myboot = function(datain,nsamples) {
  #datain: output from subatheraggfun
  observedmean = tapply(datain$Corr,INDEX = datain$MB,mean)
  n.each = c(table(datain$MB))
  observedse = c(tapply(datain$Corr,INDEX = datain$MB,sd)/sqrt(n.each))
  
  observedt = observedmean/observedse
  boott = matrix(NA,nrow=length(observedt),ncol=nsamples)
  for (i in 1:nsamples) {
    newid = sample(1:nsubject,replace=TRUE)
    stop('this code has errors')
    newdata = datain[datain$id%in%newid,]
    newmean = tapply(newdata$Corr,INDEX = newdata$MB,mean)
    newse = tapply(newdata$Corr,INDEX = newdata$MB,sd)/sqrt(n.each)
    boott[,i] = newmean/newse
  }
  setstat = apply(boott,1,sd)
 
  newdf = data.frame('MB.Factor' = names(observedmean), 'Mean.Corr' = observedmean, 'SE.Corr'=observedse,'tstat'=observedt,'SE.tstat'=setstat)
  
  newdf$MB.Factor = factor(newdf$MB.Factor,levels=c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12"))
  
  newdf
}




heatmap.pownet <- function(Image, ylab = "", xlab="",diag=F, upper=0.5,lower=-0.5)
{
  ##Read in Power Atlas
  #atlas=read.csv("~/Dropbox/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged.csv")
  
  neworder = read.csv("./OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv",header=FALSE)[,1]
  # labs=names(table(atlas$CommunityLabel)) ##Retrieve community names
  #roi=atlas$ROI[order(atlas$CommunityLabel)]       ##This will provide the order of the ROIs that places them in proper communities
  
  count.at=table(neworder)       ##How many ROIs are in each community?
  
  
  if(diag) Image= Image-diag(diag(Image))
  
  Grid.at = cumsum(count.at)
  
  #comm.labs=atlas[roi,"CommunityLabel"]
  #labs=comm.labs[Grid.at]
  labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA","Uncr")
  ticks=vector()
  for (i in 1:14){
    if (i==1) ticks[i]=Grid.at[i]/2
    else ticks[i]=(Grid.at[i]+Grid.at[i-1])/2
  }
  
  #lower = mean(Image) - 3 * sd(Image)
  #upper = mean(Image) + 3 * sd(Image)

  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

  Image[Image < lower] = lower
  Image[Image > upper] = upper
  diag(Image) = NA
  image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
        zlim = c(lower, upper), axes = FALSE, col = jet.colors(100), 
        xlab = xlab, ylab = ylab)
  
  #title(main=main)
  
  #abline(v=c(0,Grid.at,264)+.5,
  #         h=c(0,Grid.at,264)+.5)
  axis(side=2,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  #axis(side=1,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  text(x=ticks, y = par("usr")[3]-7, labels=labs2,pos=1,las=1,cex.axis=0.5,srt=45,xpd=TRUE,adj=1) 
}


## Create a plotting function that does not include nodes labeled as uncertain:
heatmap.pownet.nouncertain <- function(Image, ylab = "", xlab="",diag=F, upper=0.5,lower=-0.5)
{
  ##Read in Power Atlas
  #atlas=read.csv("~/Dropbox/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged.csv")
  
  neworder = read.csv("./OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv",header=FALSE)[,1]
  # labs=names(table(atlas$CommunityLabel)) ##Retrieve community names
  #roi=atlas$ROI[order(atlas$CommunityLabel)]       ##This will provide the order of the ROIs that places them in proper communities
  Image = Image[neworder!=14,neworder!=14]
  neworder=neworder[neworder!=14]
  count.at=table(neworder)       ##How many ROIs are in each community?
  
  
  if(diag) Image= Image-diag(diag(Image))
  
  Grid.at = cumsum(count.at)
  
  #comm.labs=atlas[roi,"CommunityLabel"]
  #labs=comm.labs[Grid.at]
  #labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA","Uncr")
  labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA")
  ticks=vector()
  for (i in 1:13){
    if (i==1) ticks[i]=Grid.at[i]/2
    else ticks[i]=(Grid.at[i]+Grid.at[i-1])/2
  }
  
  #lower = mean(Image) - 3 * sd(Image)
  #upper = mean(Image) + 3 * sd(Image)
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  diag(Image) = NA
  image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
             zlim = c(lower, upper), axes = FALSE, col = jet.colors(100), 
             xlab = xlab, ylab = ylab)
  
  #title(main=main)
  
  #abline(v=c(0,Grid.at,264)+.5,
  #         h=c(0,Grid.at,264)+.5)
  axis(side=2,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  #axis(side=1,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  text(x=ticks, y = par("usr")[3]-7, labels=labs2,pos=1,las=1,cex.axis=0.5,srt=45,xpd=TRUE,adj=1) 
}

## Create a plotting function that does not include nodes labeled as uncertain:
heatmap.pownet.subcortcer <- function(Image, ylab = "", xlab="",diag=F, upper=0.5,lower=-0.5)
{
  ##Read in Power Atlas
  #atlas=read.csv("~/Dropbox/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged.csv")
  
  neworder = read.csv("./OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv",header = FALSE)[,1]
  # labs=names(table(atlas$CommunityLabel)) ##Retrieve community names
  #roi=atlas$ROI[order(atlas$CommunityLabel)]       ##This will provide the order of the ROIs that places them in proper communities
  Image = Image[1:17,1:17]
  neworder=neworder[1:17]
  count.at=table(neworder)       ##How many ROIs are in each community?
  
  
  if(diag) Image= Image-diag(diag(Image))
  
  Grid.at = cumsum(count.at)
  
  #comm.labs=atlas[roi,"CommunityLabel"]
  #labs=comm.labs[Grid.at]
  #labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA","Uncr")
  labs2=c("Subc","Cer")
  ticks=vector()
  for (i in 1:2){
    if (i==1) ticks[i]=Grid.at[i]/2
    else ticks[i]=(Grid.at[i]+Grid.at[i-1])/2
  }
  
  #lower = mean(Image) - 3 * sd(Image)
  #upper = mean(Image) + 3 * sd(Image)
  
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  diag(Image) = NA
  image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
             zlim = c(lower, upper), axes = FALSE, col = jet.colors(100), 
             xlab = xlab, ylab = ylab)
  
  #title(main=main)
  
  #abline(v=c(0,Grid.at,264)+.5,
  #         h=c(0,Grid.at,264)+.5)
  axis(side=2,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  #axis(side=1,at=ticks,labels=labs2,las=1,cex.axis=1,tick=F) 
  text(x=ticks, y = par("usr")[3]-1, labels=labs2,pos=1,las=1,cex.axis=0.5,srt=45,xpd=TRUE,adj=1) 
}

# df_Img <- melt(Image)
# 
# ggplot(data = df_Img,
#        aes(x = Var1, y = Var2, fill = value)) +
#   geom_raster() +
#   scale_y_discrete(name="",breaks=c(0,Grid.at)+.5,labels=labs) +
#   xlab("") +
#   ylab("") +
#   geom_vline(xintercept = c(0,Grid.at,264)+.5) +
#   geom_hline(yintercept = c(0,Grid.at,264)+.5) +
#   scale_fill_distiller(palette = "RdBu",guide = F) +
#   theme(axis.text.y = element_text(angle=45, hjust = 1),
#         axis.text.x = element_blank())
  
  #ggtitle("ggplot heatmap")
    
    # abline(v=c(0,Grid.at,264)+.5,
    #        h=c(0,Grid.at,264)+.5)



# heatmap for subcortical and cerebellum:

heatmap.generic <- function(Image, ylab = "", xlab="",diag=F, upper=0.5,lower=-0.5) {
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  
  Image[Image < lower] = lower
  Image[Image > upper] = upper
  diag(Image) = NA
  image.plot(x = 1:dim(Image)[2], y = 1:dim(Image)[1], z = t(Image), 
        zlim = c(lower, upper), axes = FALSE, col = jet.colors(100), 
        xlab = xlab, ylab = ylab)
  }



# powcircplot = function(diff.mat,lower=-0.2,upper=0.2) {
#   require(circlize)
#   #install.packages("BiocManager")
#   #BiocManager::install("ComplexHeatmap")
#   require(ComplexHeatmap)
#   
#   label = read.csv('~/Dropbox/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort_withroi.csv',header = FALSE)
#   labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA","Uncr")
#   #labs3 = c(-1,0,1,2,3,4,5,6,7,8,9,11,12,14) # codes in label
#   
#   
#   keep = which(label$V2 %in% c(-1,0,5,9,11))
#   
#   label0 = label$V2[keep]
#   #label1 = sort(label0)
#   label1 = label0
#   lab.idx = label$V1[keep]
#   #lab.idx = sort(label0,index.return=T)$ix
#   cols = c(rep("brown",sum(label1==-1)),
#            rep("cadetblue2",sum(label1==0)),
#            rep("tomato",sum(label1==5)),
#            rep("royalblue3",sum(label1==9)),
#            rep("yellow",sum(label1==11)))
#   mat = as.matrix(diff.mat[keep,keep])
#   p = nrow(mat)
#   
#   node.names = paste0("power",lab.idx)
#   node.names = factor(node.names,levels=unique(node.names))
#   
#   colormap = colorRamp2(c(lower,0,upper),c("blue","white","red"))
#   
#   circos.par("track.height" = 0.2,start.degree=90)
#   circos.initialize(factors = node.names,xlim=c(0,1))
#   circos.track(bg.col=cols,ylim=c(0,1),
#                panel.fun=function(x,y){
#                  circos.text(0.5,CELL_META$ylim[1]+uy(6,'mm'),CELL_META$sector.index,
#                              facing='clockwise',niceFacing = T,
#                              cex=0.7)
#                })
#   for(i in 2:p){
#     for(j in 1:(i-1)){
#       if(mat[i,j]!=0){
#         circos.link(node.names[j],0,node.names[i],0,col=colormap(mat[i,j]),lwd=1.5)
#       }
#     }
#   }
#   
#   lgd_links = Legend(at = c(lower,0,upper), col_fun = colormap, 
#                      title_position = "topleft", title = "Links")
#   lgd_blocks = Legend(labels=c("Subcortical","Cerebellum","Default Mode","Visual","Fronto-Parietal Task"), 
#                       # type="grid",grid_height = unit(4, "mm"), grid_width = unit(4, "mm"),
#                       legend_gp=gpar(fill=c("brown","cadetblue2","tomato","royalblue3","yellow")),title_position = "topleft",title="Modules")
#   
#   lgd_list = packLegend(lgd_links, lgd_blocks)
#   draw(lgd_list, x = unit(10, "mm"), y = unit(10, "mm"), just = c("left", "bottom"))
#   
#   circos.clear()
# }
# 



# This function takes the input array as the data
# with an eye towards eventually calculating boot se
calc.prop.active = function(z.corr.array,p.thresh,indexvariable,colLabels=NULL) {
  # assumes n = 32
  # z.corr.array is 264x264x9x32. These are Fisher z-transformed correlations
  # uses sample size correction so power is standardized across acquisitions
  # indexvariable is a vector of JPower communities corresponding to the rows
  # and columns of z.corr.array
  
  mean.z.cor = apply(z.corr.array,c(1,2,3),mean,na.rm=TRUE)
  mean.z.cor.sd = apply(z.corr.array,c(1,2,3),sd,na.rm=TRUE)
  cd.z.cor = mean.z.cor / mean.z.cor.sd
  tstat.z.cor = sqrt(32)*cd.z.cor
  
  # find the Bonferonni level:
  # use two-sided pvalue
  tThresh = abs(qt(p.thresh/2,df = 31))
  
  propSigni=apply(abs(tstat.z.cor)>tThresh,MARGIN = c(1,3), FUN=mean,na.rm=TRUE)
  # take the mean by community:
  
  nComm = length(table(indexvariable))
  propSigni_comm = matrix(0,nComm,9)
  
  
  for (k in 1:9) {
    temp = tapply(X = propSigni[,k],INDEX=indexvariable,FUN=mean,na.rm=TRUE)
    propSigni_comm[,k] = temp
  }
  tdata = data.frame('Community'=names(temp),propSigni_comm)
  names(tdata)[2:10] = colLabels

  # take the overall mean:
  allrow = apply(propSigni,2,mean,na.rm=TRUE)
  
  alldat = data.frame('Community'='All',t(allrow))
  names(alldat)[2:10] = colLabels
  rbind(tdata,alldat)
}


calc.cohensd.comm = function(z.corr.array,indexvariable,colLabels=NULL) {
  # assumes n = 32
  # z.corr.array is 264x264x9x32. These are Fisher z-transformed correlations
  # uses sample size correction so power is standardized across acquisitions
  # indexvariable is a vector of JPower communities corresponding to the rows
  # and columns of z.corr.array
  
  mean.z.cor = apply(z.corr.array,c(1,2,3),mean,na.rm=TRUE)
  mean.z.cor.sd = apply(z.corr.array,c(1,2,3),sd,na.rm=TRUE)
  cd.z.cor = mean.z.cor / mean.z.cor.sd
  
 
  mean.cd=apply(abs(cd.z.cor),MARGIN = c(1,3), FUN=mean, na.rm=TRUE)

  nComm = length(table(indexvariable))
  cd_comm = matrix(0,nComm,9)
  
  # here, take the mean cohen's d for a node in the community
  for (k in 1:9) {
    temp = tapply(X = mean.cd[,k],INDEX=indexvariable,FUN=mean)
    cd_comm[,k] = temp
  }
  tdata = data.frame('Community'=names(temp),cd_comm)
  names(tdata)[2:10] = colLabels
  
  # take the overall mean:
  allrow = apply(mean.cd,2,mean,na.rm=TRUE)
  
  alldat = data.frame('Community'='All',t(allrow))
  names(alldat)[2:10] = colLabels
  rbind(tdata,alldat)
}

#####################################
create.data.to.compare = function(alldata,subject_exclude_df,acquisition1,acquisition2) {
  # alldata: nnode x nnode x 9 x nsubject
  # acquisitions:
  #   1 sb 3.3
  #   2 sb 2
  #   4 mb 3
  #   5 mb 4
  #   6 mb 6
  #   7 mb 8
  #   8 mb 9
  #   9 mb 12
  
  # e.g., set acquisition2 to 8 for creating datasets to test whether acquisition1 differs from mb 8
  # here, we will only include a subject if they have usabel scans for both acquisitions:
  subset = alldata[,,c(acquisition1,acquisition2),]
  subset_subj_include = apply(subject_exclude_df[,c(acquisition1+1,acquisition2+1)],1,max)==0
  subset = subset[,,,subset_subj_include]
  nnode = dim(subset)[1]
  for (j in 1:nnode) {
    subset[j,j,,]=0
  }
  
  # next, if edge is missing data in acquisition1, make missing in acq2, and vice versa:
  na_acquisition = is.na(subset)
  #any_na_acquisition = apply(na_acquisition,c(1,2,4),max)
  any_na_acquisition = (na_acquisition[,,1,]+na_acquisition[,,2,])>0
  subset1 = subset[,,1,]
  subset2 = subset[,,2,]
  subset1[any_na_acquisition] = NA
  subset2[any_na_acquisition] = NA
  rm(subset)
  # check NAs match:
  #all(is.na(subset1)==is.na(subset2))
  
  nsubject = dim(subset1)[3]
  nedge = choose(nnode,2)
  array(c(subset1,subset2),dim = c(nnode,nnode,nsubject,2))
}

calc.diff.prop.activated = function(datain,p.threshold,indexvariable) {
  # datain: output from create.two.data, nnode x nnode x subject x 2
  # p.threshold: user-specified, e.g., sets to 0.05/choose(nnode,2)
  nsubject = dim(datain)[3]
  mean.z.cor = apply(datain,c(1,2,4),mean,na.rm=TRUE)
  mean.z.cor.sd = apply(datain,c(1,2,4),sd,na.rm=TRUE)
  cd.z.cor = mean.z.cor / mean.z.cor.sd
  tstat.z.cor = sqrt(nsubject)*cd.z.cor
  
  # find the Bonferonni level:
  # use the two-sided p-value
  tThresh = abs(qt(p.threshold/2,df = (nsubject-1)))
  
  propSigni=apply(abs(tstat.z.cor)>tThresh,MARGIN = c(1,3), FUN=mean,na.rm=TRUE)
  # take the mean by community:
  nComm = length(table(indexvariable))
  propSigni_comm = matrix(0,nComm,2)
  
  for (k in 1:2) {
    temp = tapply(X = propSigni[,k],INDEX=indexvariable,FUN=mean,na.rm=TRUE)
    propSigni_comm[,k] = temp
  }
  diffcomm = propSigni_comm[,1] - propSigni_comm[,2]
  
  # take the overall mean:
  allrow = apply(propSigni,2,mean,na.rm=TRUE)
  
  diffall = allrow[1]-allrow[2]
  diffout = c(diffcomm,diffall)
  names(diffout) = c(names(temp),'All')
  diffout
}


perm.calc.diff.prop.activated = function(datain,p.threshold,indexvariable,seed=NULL) {
  # datain: output from create.two.data, nnode x nnode x subject x 2
  # p.threshold: user-specified, e.g., sets to 0.05/choose(nnode,2)
  
  if(!is.null(seed)) set.seed(seed)
  nsubject = dim(datain)[3]
  perm.datain=array(0,dim(datain))
  for (i in 1:nsubject) {
    myperm=sample(1:2)
    perm.datain[,,i,] = datain[,,i,myperm]
  }
  calc.diff.prop.activated(perm.datain,p.threshold,indexvariable)  
}

