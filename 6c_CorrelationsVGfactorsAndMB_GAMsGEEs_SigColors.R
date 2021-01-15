#####################################################################
#
# This file analyzes impacts of g-factor and multiband factor for
# a subset of intra-commmunity edges (6 regions/communities)
# 1. GAMMs analyzing impact of gfactor on correlations
# 2. GEEs analyzing subset of edges for different communities and MB factor
#
######################################################################


library(R.matlab)
library(lme4) 
library(lmerTest)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(matlab)
library(mgcv)
library(itsadug)
library(geepack)
library(xtable)

setwd('~/risk_share/')
      
source('./OptimalSMS_rsfMRI/Programs/Functions/fconn_fun.R')
alldat = readMat('./OptimalSMS_rsfMRI/Results/CorrSDpower264_9p.mat')
alldat_tf = readMat('./OptimalSMS_rsfMRI/Results/CorrSDpower264_9p_tf.mat')


demog = read.csv('./OptimalSMS_rsfMRI/Data2/Subject_Demog.csv')

# subject list from matlab processing: 2_ConsolidateCorrMat_cluster.m
# Data must be in the same order used to create correlation arrays
excludeScans = read.csv('./OptimalSMS_rsfMRI/Results/exclude_subjects.csv');
subjectList=excludeScans[,1]
subjectData = data.frame('ID'=subjectList,'order' = 1:length(subjectList)) # order is the ordering in the matlab
# CorrSdpower264_v2.mat

subjectData2 = merge(subjectData,demog)
head(subjectData2)
nrow(subjectData2)

# Data must be in the order of subjectList:
subjectData2 = subjectData2[order(subjectData2$order),]
head(subjectData2)
nrow(subjectData2)

dim(alldat$allmb.cor.srt)

nsubject = dim(alldat$allmb.cor.srt)[4]
node1 = c(1:264)%*%t(rep(1,264))
node2 = rep(1,264)%*%t(c(1:264))

varnames = c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12")
varnames_zcor = paste0('zcor_',varnames)
varnames_gg = paste0('gg_',varnames)


## Master datasets used for all ROIs:
subj.cor.df = NULL
for (i in 1:nsubject) {
  temp = apply(alldat$allmb.cor.srt[,,,i],3,net2vec)
  temp.df = data.frame(atanh(temp))
  names(temp.df) = varnames_zcor
  temp.df$id = i
  temp.df$gender = subjectData2$Gender[i]
  temp.df$scanner = subjectData2$CSIC_BITC[i]
  temp.df$node1 = net2vec(node1)
  temp.df$node2 = net2vec(node2)
  ptemp = temp.df%>%pivot_longer(cols=starts_with('zcor'),names_to='Acquisition',names_prefix='zcor_',values_to='zcor')
  
  temp2 = apply(alldat$allmb.gg.srt[,,,i],3,net2vec)
  temp2.df = data.frame(temp2)
  names(temp2.df)=varnames_gg
  ptemp2=temp2.df%>%pivot_longer(cols=starts_with('gg_'),names_to='Acquisition2',names_prefix='gg_',values_to='gstar_factor')
  ptemp3 = cbind(ptemp,ptemp2[,2])
  
  subj.cor.df = rbind(subj.cor.df,ptemp3)
}
head(subj.cor.df)



subj_tf.cor.df = NULL
for (i in 1:nsubject) {
  temp = apply(alldat_tf$allmb.tf.cor.srt[,,,i],3,net2vec)
  temp.df = data.frame(atanh(temp))
  names(temp.df) = varnames_zcor
  temp.df$id = i
  temp.df$gender = subjectData2$Gender[i]
  temp.df$scanner = subjectData2$CSIC_BITC[i]
  temp.df$node1 = net2vec(node1)
  temp.df$node2 = net2vec(node2)
  ptemp = temp.df%>%pivot_longer(cols=starts_with('zcor'),names_to='Acquisition',names_prefix='zcor_',values_to='zcor')
  
  temp2 = apply(alldat_tf$allmb.tf.gg.srt[,,,i],3,net2vec)
  temp2.df = data.frame(temp2)
  names(temp2.df)=varnames_gg
  ptemp2=temp2.df%>%pivot_longer(cols=starts_with('gg_'),names_to='Acquisition2',names_prefix='gg_',values_to='gstar_factor')
  ptemp3 = cbind(ptemp,ptemp2[,2])
  
  subj_tf.cor.df = rbind(subj_tf.cor.df,ptemp3)
}
head(subj_tf.cor.df)


# Merge with original power ROI labels:
power_resort = read.csv('./OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv',header = FALSE)
names(power_resort) = c('Community','ROI','BRisk_ROI')

temp1 = data.frame('Community1'=power_resort$Community,'node1_power'=power_resort$ROI,'node1'=power_resort$BRisk_ROI)
temp2 = data.frame('Community2'=power_resort$Community,'node2_power'=power_resort$ROI,'node2'=power_resort$BRisk_ROI)

subj.cor.df = merge(subj.cor.df,temp1)
subj.cor.df = merge(subj.cor.df,temp2)
subj.cor.df$PowerNode1Node2 = paste0(subj.cor.df$node1_power,'_',subj.cor.df$node2_power)
subj.cor.df$Acquisition=factor(subj.cor.df$Acquisition)
subj.cor.df$Acquisition=relevel(subj.cor.df$Acquisition,ref = 'SB.2.mm')

subj_tf.cor.df = merge(subj_tf.cor.df,temp1)
subj_tf.cor.df = merge(subj_tf.cor.df,temp2)
subj_tf.cor.df$PowerNode1Node2 = paste0(subj_tf.cor.df$node1_power,'_',subj.cor.df$node2_power)
subj_tf.cor.df$Acquisition=factor(subj_tf.cor.df$Acquisition)
subj_tf.cor.df$Acquisition=relevel(subj_tf.cor.df$Acquisition,ref = 'SB.2.mm')



# These are the node IDs from the original sorting: 
# Subcortical: Right Thalamus, Left Thalamus: 225, 224
# Subcortical: Right Putamen, Left Putamen: 230 227
# Cerebellum: Right Cerebellum, Left Cerebellum: 245 244
# Executive control: 196, 201, 192, 195 (in HBM paper)
# Salience: 217, 212, 220, 218
# DMN: 110, 92, 88, 109

# selected nodes:
selnode1 = c(225,230,245,201,201,201,196,196,195,220,220,220,218,218,217,110,110,110,109,109,92)
selnode2 = c(224,227,244,196,195,192,195,192,192,218,217,212,217,212,212,109,92,88,92,88,88)
sel_PowerNode1Node2=paste0(selnode1,'_',selnode2)

subset.cor.df = subj.cor.df[subj.cor.df$PowerNode1Node2%in%sel_PowerNode1Node2,]
# Check whether we grabbed all the edges we wanted. Should be 21 edges:
length(table(subset.cor.df$PowerNode1Node2)) # Correctly subset.

subset_tf.cor.df = subj_tf.cor.df[subj_tf.cor.df$PowerNode1Node2%in%sel_PowerNode1Node2,]
length(table(subset_tf.cor.df$PowerNode1Node2))

# Figure that illustrate correlations get squished with higher g-factor:
subset.cor.df%>%ggplot(aes(x=gstar_factor,y=zcor))+geom_point()


subset.cor.df$PowerNode1Node2=factor(subset.cor.df$PowerNode1Node2)
contrasts(subset.cor.df$PowerNode1Node2)=contr.sum


subset.cor.df$id_factor = factor(subset.cor.df$id)
# interaction between subject id and power node:
subset.cor.df$id_PowerNode1Node2 = factor(paste0(subset.cor.df$id,'_',subset.cor.df$PowerNode1Node2))

subset_tf.cor.df$id_factor=factor(subset_tf.cor.df$id)
subset_tf.cor.df$id_PowerNode1Node2 = factor(paste0(subset_tf.cor.df$id,'_',subset_tf.cor.df$PowerNode1Node2))


# full model:
ggfactor_model = gam(zcor~s(gstar_factor,k=-1)+s(id_factor,bs='re')+s(id_PowerNode1Node2,bs='re')+scanner+gender+PowerNode1Node2,data=subset.cor.df)
summary(ggfactor_model)
plot(ggfactor_model)
hist(residuals(ggfactor_model))

# plot_smooth requires you to choose a level for each factor. Choose central node: 
median(tapply(subset.cor.df$zcor,subset.cor.df$PowerNode1Node2,FUN=mean,na.rm=TRUE))
# use 109_92


pdf(file='~/risk_share/OptimalSMS_rsfMRI/Documents/Figures/ZCor_Gfactor_gamm.pdf')
plot_smooth(ggfactor_model,view='gstar_factor',rm.ranef=TRUE, rug=FALSE, xlab='g*-factor',ylab='Z Correlation',cex.lab=1.3,xlim=c(0.8,8),cond=list(PowerNode1Node2='109_92'),ylim=c(-0.1,0.4))
dev.off()

# calculate Cohen's D for each edge:
grouped = subset.cor.df%>%group_by(PowerNode1Node2,Acquisition)
subset.cd.df = summarise(grouped, mean_zcor=mean(zcor,na.rm=TRUE), sd_zcor=sd(zcor,na.rm=TRUE), mean_gg=mean(gstar_factor,na.rm=TRUE))
subset.cd.df$cd = subset.cd.df$mean_zcor/subset.cd.df$sd_zcor


# This result is currently not in manuscript. Would have to construct a bootstrap test of significance.
# These confidence intervals are not correct:
plot(cd~mean_gg,data=subset.cd.df)
ggfactor_cd_model = gam(cd~s(mean_gg,k=-1)+PowerNode1Node2,data=subset.cd.df)
plot_smooth(ggfactor_cd_model,view='mean_gg',rug=FALSE, xlab='gg factor',ylab='Cohen\'s D',cex.lab=1.3,xlim=c(0.8,10))
# negative trend 
######################

  
  
#############################
# Intra-community correlations
########################

# Subcortical: Right Thalamus, Left Thalamus: 225, 224
# Subcortical: Right Putamen, Left Putamen: 230 227
# Cerebellum: Right Cerebellum, Left Cerebellum: 245 244
# Executive control: 196, 201, 192, 195 (in HBM paper)
# Salience: 217, 212, 220, 218
# DMN: 110, 92, 88, 109
subset.cor.df = subset.cor.df[order(subset.cor.df$id),]
subset_tf.cor.df = subset_tf.cor.df[order(subset_tf.cor.df$id),]
  
subset.cor.df$Acquisition = relevel(subset.cor.df$Acquisition,ref="SB.2.mm")
subset_tf.cor.df$Acquisition = relevel(subset_tf.cor.df$Acquisition,ref="SB.2.mm")
  
thalamus =  subset.cor.df[subset.cor.df$PowerNode1Node2=='225_224',]
thalamus_tf = subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='225_224',]

putamen = subset.cor.df[subset.cor.df$PowerNode1Node2=='230_227',]
putamen_tf=subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='230_227',]

cerebellum=subset.cor.df[subset.cor.df$PowerNode1Node2=='245_244',]
cerebellum_tf=subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='245_244',]

executive=subset.cor.df[subset.cor.df$node1_power%in%c(192,195,196,201) & subset.cor.df$node2_power%in%c(192,195,196,201),]
executive_tf = subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(192,195,196,201) & subset_tf.cor.df$node2_power%in%c(192,195,196,201),]

salience=subset.cor.df[subset.cor.df$node1_power%in%c(212,217,218,220) & subset.cor.df$node2_power%in%c(212,217,218,220),]
salience_tf=subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(212,217,218,220) & subset_tf.cor.df$node2_power%in%c(212,217,218,220),]

dmn = subset.cor.df[subset.cor.df$node1_power%in%c(110,92,88,109) & subset.cor.df$node2_power%in%c(110,92,88,109),]  
dmn_tf = subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(110,92,88,109) & subset_tf.cor.df$node2_power%in%c(110,92,88,109),]

###
# thalamus.model = lmer(zcor~Acquisition+gender+scanner+(1|id),data=thalamus)  
# thalamus_tf.model = lmer(zcor~Acquisition+gender+scanner+(1|id),data=thalamus_tf) 
# 
# putamen.model = lmer(zcor~Acquisition+gender+scanner+(1|id),data=putamen)
# putamen_tf.model=lmer(zcor~Acquisition+gender+scanner+(1|id),data=putamen_tf)
# 
# cerebellum.model=lmer(zcor~Acquisition+gender+scanner+(1|id),data=cerebellum)
# cerebellum_tf.model=lmer(zcor~Acquisition+gender+scanner+(1|id),data=cerebellum_tf)
# 
# executive.model = lmer(zcor~Acquisition+gender+scanner+PowerNode1Node2+(1|id/PowerNode1Node2),data=executive)
# executive_tf.model = lmer(zcor~Acquisition+gender+scanner+PowerNode1Node2+(1|id/PowerNode1Node2),data=executive_tf)
# 
# salience.model = lmer(zcor~Acquisition+gender+scanner+(1|id/PowerNode1Node2),data=salience)
# salience_tf.model = lmer(zcor~Acquisition+gender+scanner+(1|id/PowerNode1Node2),data=salience_tf)
# 
# dmn.model = lmer(zcor~Acquisition+gender+scanner+(1|id/PowerNode1Node2),data=dmn)
# dmn_tf.model = lmer(zcor~Acquisition+gender+scanner+(1|id/PowerNode1Node2),data=dmn_tf)

for.table.lmer = function(out.model,region) {
  temp = summary(out.model)$coefficients[1:9,c(1,2,5)]
  temp2 = cbind(round(temp[,1:2],2),round(temp[,3],3))
  temp3 = c(paste0(temp2[1,1],'(',temp2[1,2],')'),paste0(temp2[2:9,1],'(',temp2[2:9,2],') ',ifelse(temp2[2:9,3]==0,'<0.001',temp2[2:9,3])))
  temp.data = data.frame(temp3)
  names(temp.data)=region
  row.names(temp.data)=c('Int.','SB 3.3mm','MB 2','MB 3','MB 4','MB 6','MB 8','MB 9','MB 12')
  temp.data
}

for.table.geepack = function(out.model,region) {
  temp = summary(out.model)$coefficients[c(1,9,3:8,2),c(1,2,4)]
  temp2 = cbind(round(temp[,1:2],2),round(temp[,3],3))
  temp3 = c(paste0(temp2[1,1],'(',temp2[1,2],')'),paste0(temp2[2:9,1],'(',temp2[2:9,2],') ',ifelse(temp2[2:9,3]==0,'$<$.001',temp2[2:9,3])))
  temp.data1 = data.frame(temp3)
  temp.data = temp.data1
  temp.data[,1] = as.character(temp.data1[,1])
  for(i in 2:9){
    if(temp2[i,3]<0.01 & sign(temp2[i,1])==1){
      temp.data[i,1] = paste0('\\textcolor{blue}{\\bf{',temp.data1[i,1], '}}')
    } else if (temp2[i,3]<0.01 & sign(temp2[i,1])==-1){
      temp.data[i,1] = paste0('\\textcolor{red}{\\bf{',temp.data1[i,1], '}}')
    }
  } 
  names(temp.data)=region
  row.names(temp.data)=c('Int.','SB3.3','MB 2','MB 3','MB 4','MB 6','MB 8','MB 9','MB 12')
  temp.data
}

#lmm.results.9p = cbind(for.table(thalamus.model,'Thalamus - 9p'),for.table(putamen.model,'Putamen - 9p'),for.table(cerebellum.model,'Cerebellum - 9p'),for.table(executive.model,'Executive - 9p'),for.table(salience.model,'Salience - 9p'),for.table(dmn.model,'DMN - 9p'))
#print(xtable(lmm.results.9p))
#lmm.results.9p.tf = cbind(for.table(thalamus_tf.model,'Thalamus - 9p+bp'),for.table(putamen_tf.model,'Putamen - 9p+bp'),for.table(cerebellum_tf.model,'Cerebellum - 9p+bp'),for.table(executive_tf.model,'Executive - 9p+bp'),for.table(salience_tf.model,'Salience - 9p+bp'),for.table(dmn_tf.model,'DMN - 9p+bp'))
#print(xtable(lmm.results.9p.tf))


################################
# switched to GEEs to have robust se (due to possible heteroscedasticity). Does not appear to impact results:
###
thalamus.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=thalamus)  
thalamus_tf.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=thalamus_tf) 

putamen.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=putamen)
putamen_tf.model=geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=putamen_tf)

cerebellum.model=geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=cerebellum)
cerebellum_tf.model=geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=cerebellum_tf)

executive.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=executive)
executive_tf.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=executive_tf)

salience.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=salience)
salience_tf.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=salience_tf)

dmn.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=dmn)
dmn_tf.model = geeglm(zcor~Acquisition+gender+scanner,id=id,corstr='exchangeable',data=dmn_tf)


gee.results.9p = cbind(for.table.geepack(thalamus.model,'Thalamus - 9p'),for.table.geepack(putamen.model,'Putamen - 9p'),for.table.geepack(cerebellum.model,'Cerebellum - 9p'),for.table.geepack(executive.model,'Executive - 9p'),for.table.geepack(salience.model,'Salience - 9p'),for.table.geepack(dmn.model,'DMN - 9p'))
print(xtable(gee.results.9p),type='latex',sanitize.text.function = function(x) {x})


gee.results.9p.tf = cbind(for.table.geepack(thalamus_tf.model,'Thalamus - 9p+bp'),for.table.geepack(putamen_tf.model,'Putamen - 9p+bp'),for.table.geepack(cerebellum_tf.model,'Cerebellum - 9p+bp'),for.table.geepack(executive_tf.model,'Executive - 9p+bp'),for.table.geepack(salience_tf.model,'Salience - 9p+bp'),for.table.geepack(dmn_tf.model,'DMN - 9p+bp'))
print(xtable(gee.results.9p.tf),type='latex',sanitize.text.function = function(x) {x})




## Re-create datasets to have desired ordering for plotting:
subset.cor.df$Acquisition = factor(subset.cor.df$Acquisition,levels=c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12"))
subset_tf.cor.df$Acquisition = factor(subset_tf.cor.df$Acquisition,levels=c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12"))

#subset.cor.df$Acquisition[subset.cor.df$Acquisition=='SB.3.3.mm']="SB3.3mm"
#subset.cor.df$Acquisition[subset.cor.df$Acquisition=='SB.2.mm']="SB2mm"
#subset_tf.cor.df$Acquisition[subset_tf.cor.df$Acquisition=='SB.3.3.mm']="SB.3.3mm"
#subset_tf.cor.df$Acquisition[subset_tf.cor.df$Acquisition=='SB.2.mm']="SB.2mm"


thalamus =  subset.cor.df[subset.cor.df$PowerNode1Node2=='225_224',]
thalamus_tf = subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='225_224',]
putamen = subset.cor.df[subset.cor.df$PowerNode1Node2=='230_227',]
putamen_tf=subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='230_227',]
cerebellum=subset.cor.df[subset.cor.df$PowerNode1Node2=='245_244',]
cerebellum_tf=subset_tf.cor.df[subset_tf.cor.df$PowerNode1Node2=='245_244',]
executive=subset.cor.df[subset.cor.df$node1_power%in%c(192,195,196,201) & subset.cor.df$node2_power%in%c(192,195,196,201),]
executive_tf = subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(192,195,196,201) & subset_tf.cor.df$node2_power%in%c(192,195,196,201),]
salience=subset.cor.df[subset.cor.df$node1_power%in%c(212,217,218,220) & subset.cor.df$node2_power%in%c(212,217,218,220),]
salience_tf=subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(212,217,218,220) & subset_tf.cor.df$node2_power%in%c(212,217,218,220),]
dmn = subset.cor.df[subset.cor.df$node1_power%in%c(110,92,88,109) & subset.cor.df$node2_power%in%c(110,92,88,109),]  
dmn_tf = subset_tf.cor.df[subset_tf.cor.df$node1_power%in%c(110,92,88,109) & subset_tf.cor.df$node2_power%in%c(110,92,88,109),]

######Settings for ggplot figure
theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)

pal <- wes_palette("Zissou1", 100, type = "continuous")
# make max gg_MB.12 = 9 to facilitate plotting


p0 = ggplot(thalamus, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Thalamus - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")


p0_tf = ggplot(thalamus_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Thalamus - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")


p1 = ggplot(putamen, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Putamen - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")

p1_tf = ggplot(putamen_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Putamen - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")

p2 = ggplot(cerebellum, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Cerebellum - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")
p2_tf = ggplot(cerebellum_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),aes(color=gstar_factor))+ggtitle('Cerebellum - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")

p3 = ggplot(executive, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('Executive - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")
p3_tf = ggplot(executive_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('Executive - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")

p4 = ggplot(salience, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('Salience - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")
p4_tf = ggplot(salience_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('Salience - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")

p5 = ggplot(dmn, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('DMN - 9p')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")
p5_tf = ggplot(dmn_tf, aes(x=Acquisition, y=zcor))+geom_boxplot(outlier.shape=NA) + geom_jitter(shape=16,position=position_jitter(0.2),size=0.75,aes(color=gstar_factor))+ggtitle('DMN - 9p+bandpass')+ylab('Z Correlation')+scale_x_discrete(labels=c("SB.3.3.mm" = "SB3.3mm", "SB.2.mm" = "SB2mm","MB.2" = "MB2", "MB.3" = "MB3", "MB.4" = "MB4", "MB.6"="MB6","MB.8" = "MB8", "MB.9" = "MB9", "MB.12" = "MB12"))+scale_color_gradientn(colours=pal,limits=c(0.9,6),oob = scales::squish,name = "g*")



png(file='./OptimalSMS_rsfMRI/Documents/Figures/Boxplots_zcorrelations_9p.png',width=13,height=12,units='in',res=500)
grid.arrange(p0,p1,p2,p3,p4,p5,nrow=3)
dev.off()

png(file='./OptimalSMS_rsfMRI/Documents/Figures/Boxplots_zcorrelations_9p_tf.png',width=13,height=12,units='in',res=500)
grid.arrange(p0_tf,p1_tf,p2_tf,p3_tf,p4_tf,p5_tf,nrow=3)
dev.off()
