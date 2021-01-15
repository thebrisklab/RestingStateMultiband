#############################################
#### Benjamin Risk
#### Edited by Raphiel Murden
#### 1. Create plots of correlation matrices
#### 2. Create plots of MB correlation ~ SB correlation
####       colored by g*-factor
#### 3. Proportion of edges that are active (barplots)
#### 4. Mean absolute Cohen's d (barplots)
##############################################

library(R.matlab)
library(lme4) 
library(lmerTest)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(matlab)
library(fields)
library(tidyr)
library(wesanderson)
library(ggpubr)


# Change working directory depending on whether running on cluster or locally:
setwd('~/risk_share')

source('./OptimalSMS_rsfMRI/Programs/Functions/fconn_fun.R')
alldat = readMat('./OptimalSMS_rsfMRI/Results/CorrSDpower264_9p.mat')
alldat_tf = readMat('./OptimalSMS_rsfMRI/Results/CorrSDpower264_9p_tf.mat')

#alldat_aroma = readMat('./OptimalSMS_rsfMRI/Results/CorrSDpower264_9p_aroma.mat')
# Note on aroma: AROMA remove a lot of signal, resulting in correlation matrices that 
#   have much reduced correlations. This is an important avenue for future research,
#   but is too much to address in current manuscript. 



######################
######################
######################
# Communities
# note: Use node2 because took lower triangle
# use 1:13 for subcortical
# use 1:17 for subcortical + cerebellum
# use 18:47 for sensory/motor
# use 80:137 for DMN
subcortical = 1:13
cerebellum = 14:17
sensomotorH = 18:47
sensomotorM = 48:52
cingulo = 53:66
auditory = 67:79
dmn = 80:137
memory = 138:142
visual = 143:173
fptaskcontr=174:198
salience=199:216
ventralatt=217:225
dorsatt=226:236
uncertain=237:264
labs2=c("Subc","Cer","SomH","SomM","CO-C","Aud","DMN","Mem","Vis","FP-C","Sal","VenA","DorA","Uncr")

# correspond to sorted indices:
labelsLong=c(rep("Subc",13),rep("Cer",4),rep("SomH",30),rep("SomM",5),rep("CO-TC",14),rep('Auditory',13),rep('Default Mode',58),rep('Memory',5),rep('Visual',31),rep('FP-TC',15),rep('Salience',18),rep("VentralAtt",9),rep("DorsalAtt",11),rep("Uncertain",38))
  
varnames = c("SB.3.3.mm","SB.2.mm", "MB.2", "MB.3","MB.4","MB.6","MB.8","MB.9","MB.12")
varnamesplots = c("1-3.3","1-2", "2", "3","4","6","8","9","12")
varnamesplots2 = c("SB 3.3mm","SB 2mm", "MB 2", "MB 3","MB 4","MB 6","MB 8","MB 9","MB 12")

################
###############

dim(alldat$allmb.cor.srt)

nsubject = dim(alldat$allmb.cor.srt)[4]
node1 = c(1:264)%*%t(rep(1,264))
node2 = rep(1,264)%*%t(c(1:264))


allmb.cor.srt = alldat$allmb.cor.srt  ##Contains communities sorted so that subcortical is first and uncertain last
allmb.gg.srt = alldat$allmb.gg.srt
allmb.sd.srt = alldat$allmb.sdsd.srt

allmb_tf.cor.srt = alldat_tf$allmb.tf.cor.srt
allmb_tf.gg.srt = alldat_tf$allmb.tf.gg.srt
allmb_tf.sd.srt = alldat_tf$allmb.tf.sdsd.srt

# allmb_aroma.cor.srt = alldat_aroma$allmb.aroma.cor.srt
# allmb_aroma.gg.srt = alldat_aroma$allmb.aroma.gg.srt
# allmb_aroma.sd.srt = alldat_aroma$allmb.aroma.sdsd.srt


# Create mean matrices for 9p:
allmb.z.cor.srt = atanh(allmb.cor.srt)
allmb.z.cor.srt[is.infinite(allmb.z.cor.srt)]=NaN
mean.z.cor = apply(allmb.z.cor.srt,c(1,2,3),mean,na.rm=TRUE)
mean.z.cor.sd = apply(allmb.z.cor.srt,c(1,2,3),sd,na.rm=TRUE)
cd.z.cor = mean.z.cor / mean.z.cor.sd

# Create mean matrices for 9p plus tf:
allmb_tf.z.cor.srt = atanh(allmb_tf.cor.srt)
allmb_tf.z.cor.srt[is.infinite(allmb_tf.z.cor.srt)]=NaN
mean_tf.z.cor = apply(allmb_tf.z.cor.srt,c(1,2,3),mean,na.rm=TRUE)
mean_tf.z.cor.sd = apply(allmb_tf.z.cor.srt,c(1,2,3),sd,na.rm=TRUE)
cd_tf.z.cor = mean_tf.z.cor / mean_tf.z.cor.sd

# # Create mean matrices for aroma:
# allmb_aroma.z.cor.srt = atanh(allmb_aroma.cor.srt)
# allmb_aroma.z.cor.srt[is.infinite(allmb_tf.z.cor.srt)]=NaN
# mean_aroma.z.cor=apply(allmb_aroma.z.cor.srt,c(1,2,3),mean,na.rm=TRUE)
# mean_aroma.z.cor.sd=apply(allmb_aroma.z.cor.srt,c(1,2,3),sd,na.rm=TRUE)
# cd_aroma.z.cor=mean_aroma.z.cor/mean_aroma.z.cor.sd

# FOR SB 3.3mm, get rid of nodes 4 and 10, corresponding 240 and 246 in resort
mean.z.cor[c(240,246),,1]=NA
mean.z.cor[c(240,246),,1]=NA
mean.z.cor[,c(240,246),1]=NA
mean.z.cor[,c(240,246),1]=NA
cd.z.cor[c(240,246),,1]=NA
cd.z.cor[c(240,246),,1]=NA
cd.z.cor[,c(240,246),1]=NA
cd.z.cor[,c(240,246),1]=NA

 
mean_tf.z.cor[c(240,246),,1]=NA
mean_tf.z.cor[c(240,246),,1]=NA
mean_tf.z.cor[,c(240,246),1]=NA
mean_tf.z.cor[,c(240,246),1]=NA
cd_tf.z.cor[c(240,246),,1]=NA
cd_tf.z.cor[c(240,246),,1]=NA
cd_tf.z.cor[,c(240,246),1]=NA
cd_tf.z.cor[,c(240,246),1]=NA

# mean_aroma.z.cor[c(240,246),,1]=NA
# mean_aroma.z.cor[c(240,246),,1]=NA
# mean_aroma.z.cor[,c(240,246),1]=NA
# mean_aroma.z.cor[,c(240,246),1]=NA
# cd_aroma.z.cor[c(240,246),,1]=NA
# cd_aroma.z.cor[c(240,246),,1]=NA
# cd_aroma.z.cor[,c(240,246),1]=NA
# cd_aroma.z.cor[,c(240,246),1]=NA


# create a single matrix with not temporally filtered above diagonal and temporally filtered below
placeholder=matrix(NA,nrow=264,ncol=264)
mean_notf_tf.z.cor = array(NA,dim=dim(mean.z.cor))
cd_notf_tf.z.cor = mean_notf_tf.z.cor
for (i in 1:9) {
  temp1=placeholder
  temp_notf=mean.z.cor[,,i]
  temp_tf=mean_tf.z.cor[,,i]
  temp1[lower.tri(temp1)]=temp_notf[lower.tri(temp1)]
  temp1[upper.tri(temp1)]=temp_tf[upper.tri(temp1)]
  mean_notf_tf.z.cor[,,i]=temp1
  
  temp_notf=cd.z.cor[,,i]
  temp_tf=cd_tf.z.cor[,,i]
  temp1[lower.tri(temp1)]=temp_notf[lower.tri(temp1)]
  temp1[upper.tri(temp1)]=temp_tf[upper.tri(temp1)]
  cd_notf_tf.z.cor[,,i]=temp1
}
# note: R plots rows with indices arranged like numbers on an axis, so that 
# in the code above, "lower.tri" actually becomes the upper triangle in plots
# This can be confirmed by commenting out the lower tri, and seeing
# that only the upper triangle appears in the plots below. 

# ############################################
# ## 9p and 9p with temporal filtering:
png(file='./OptimalSMS_rsfMRI/Documents/Figures/ZCorr_9p_notf_tf.png',width=14,height=12,units='in',res=1000)
par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
for(i in 1:9) {
  heatmap.pownet.nouncertain(mean_notf_tf.z.cor[,,i],lower=-0.6,upper = 0.6)
  title(varnamesplots2[i],cex.main=2)
}
dev.off()

png(file='./OptimalSMS_rsfMRI/Documents/Figures/ZCorr_subcort_cereb_9p_notf_tf.png',width=14,height=12,units='in',res=1000)
par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
for(i in 1:9) {
  heatmap.pownet.subcortcer(mean_notf_tf.z.cor[1:17,1:17,i],lower=-0.6,upper = 0.6)
  title(varnamesplots2[i],cex.main=2)
}
dev.off()

png(file='./OptimalSMS_rsfMRI/Documents/Figures/CohensD_9p_notf_tf.png',width=14,height=12,units='in',res=1000)
par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
for (i in 1:9) {
  heatmap.pownet.nouncertain(cd_notf_tf.z.cor[,,i],lower=-2.5,upper=2.5)
  title(varnamesplots2[i],cex.main=2)
}
dev.off()

png(file='./OptimalSMS_rsfMRI/Documents/Figures/CohensD_subcort_cereb_9p_notf_tf.png',width=14,height=12,units='in',res=1000)
par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
for (i in 1:9) {
  heatmap.pownet.subcortcer(cd_notf_tf.z.cor[1:17,1:17,i],lower=-2.5,upper=2.5)
  title(varnamesplots2[i],cex.main=4)
}
dev.off()

# Look at whether effect sizes decrease with temporal filtering:
apply(abs(cd.z.cor),3,mean,na.rm=TRUE)
apply(abs(cd_tf.z.cor),3,mean,na.rm=TRUE)
# overall there is a decrease.

# 
# ##############################
# ##############################
# # AROMA:
# 
# png(file='./OptimalSMS_rsfMRI/Documents/Figures/ZCorr_aroma_aroma.png',width=14,height=12,units='in',res=1000)
# par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
# for(i in 1:9) {
#   heatmap.pownet.nouncertain(mean_aroma.z.cor[,,i],upper = 0.6,lower=-0.6)
#   title(varnamesplots2[i],cex.main=2)
# }
# dev.off()
# 
# png(file='./OptimalSMS_rsfMRI/Documents/Figures/ZCorr_subcort_cereb_aroma_aroma.png',width=14,height=12,units='in',res=1000)
# par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
# for (i in 1:9) {
#   heatmap.pownet.subcortcer(mean_aroma.z.cor[1:17,1:17,i],lower=-0.6,upper=0.6)
#   title(varnamesplots2[i],cex.main=4)
# }
# dev.off()
# 
# png(file='./OptimalSMS_rsfMRI/Documents/Figures/CohensD_aroma_aroma.png',width=14,height=12,units='in',res=1000)
# par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
# for(i in 1:9) {
#   heatmap.pownet.nouncertain(cd_aroma.z.cor[,,i],lower=-2.5,upper=2.5)
#   title(varnamesplots2[i],cex.main=2)
# }
# dev.off()
# 
# 
# png(file='./OptimalSMS_rsfMRI/Documents/Figures/CohensD_subcort_cereb_aroma_aroma.png',width=14,height=12,units='in',res=1000)
# par(mfrow=c(3,3),mar=c(4,7,3,2),oma=c(0,0,0,1))
# for (i in 1:9) {
#   heatmap.pownet.subcortcer(cd_aroma.z.cor[1:17,1:17,i],lower=-2.5,upper=2.5)
#   title(varnamesplots2[i],cex.main=4)
# }
# dev.off()
# 
# 
# 
# 
# 
###################################
###########################
# Create scatterplots of correlations for different multiband factors

# Code below creates a dataset of vectorized correlation matrices that can be plotted using ggplot.
# This will need to be repeated for temporal filtering (mean_tf.z.cor, allmb_tf.sd.srt, etc)

####################Construct dataset with NO BANDPASS
temp = apply(mean.z.cor,3,net2vec)
netdf_9p = data.frame(temp)
names(netdf_9p) = paste0('z_',varnames)

netdf_9p$node1 = net2vec(node1)
netdf_9p$node2 = net2vec(node2)
mean.sd = apply(apply(allmb.sd.srt,c(1,2,3),mean,na.rm=TRUE),3,net2vec)
temp = data.frame(mean.sd)
names(temp) = paste0('sdsd_',varnames)
netdf_9p = cbind(netdf_9p,temp)

mean.gg = apply(apply(allmb.gg.srt,c(1,2,3),mean,na.rm=TRUE),3,net2vec)
temp=data.frame(mean.gg)
names(temp) = paste0('gg_',varnames)
netdf_9p=cbind(netdf_9p,temp)

temp = data.frame(apply(cd.z.cor,3,net2vec))
names(temp) = paste0('cd_',varnames)
netdf_9p=cbind(netdf_9p,temp)

####################Construct dataset WITH BANDPASS
temp = apply(mean_tf.z.cor,3,net2vec)
netdf_9p_tf = data.frame(temp)
names(netdf_9p_tf) = paste0('z_',varnames)

netdf_9p_tf$node1 = net2vec(node1)
netdf_9p_tf$node2 = net2vec(node2)
mean.sd = apply(apply(allmb_tf.sd.srt,c(1,2,3),mean,na.rm=TRUE),3,net2vec)
temp = data.frame(mean.sd)
names(temp) = paste0('sdsd_',varnames)
netdf_9p_tf = cbind(netdf_9p_tf,temp)

mean.gg_tf = apply(apply(allmb_tf.gg.srt,c(1,2,3),mean,na.rm=TRUE),3,net2vec)
temp=data.frame(mean.gg_tf)
names(temp) = paste0('gg_',varnames)
netdf_9p_tf=cbind(netdf_9p_tf,temp)

temp = data.frame(apply(cd_tf.z.cor,3,net2vec))
names(temp) = paste0('cd_',varnames)
netdf_9p_tf=cbind(netdf_9p_tf,temp)

save(netdf_9p,netdf_9p_tf,file = '~/risk_share/OptimalSMS_rsfMRI/Results/vecnet_dataframe.RData')




######Settings for ggplot figure
theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)


pal <- wes_palette("Zissou1", 100, type = "continuous")
# make max gg_MB.12 = 9 to facilitate plotting

quantile(netdf_9p$gg_MB.12,0.99)



# Create 8x2 plots with z_SB.2.mm on x-axis in all plots, and y-axis SB.3.3.mm, MB 2, MB3, MB4, MB6, MB8, MB9, MB12
# fix axis labels, "SB 2 mm" "MB 12" , legend label, "g*"
# tips on multiple plots I have found helpful
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

# the first column will be the 8 plots for the 9p processing (no temporal filtering)
# the second column will be the 8 plots for the 9p+temporal filtering (i.e., bandpass)
# Create another set of figures for the bandpass filtered data. In these bandpass filtered data, we will used the new gg_MB from 
# the temporally filtered data. Use the same limits (4), in which I expect all points to be pretty solidly blue. 
# NOTE: some outliers in SB.3.3, define axis limits to avoid these (they aren't analyzed, it's a registration thing)

# note: i am able to use all points on the cluster
#mysample=sample(1:nrow(netdf_9p),10000)
#netdf_9p_sample = netdf_9p[mysample,]

xy.ratio = 3/4

sp_33 = ggplot(netdf_9p, aes(z_SB.2.mm,z_SB.3.3.mm))+
  geom_point(aes(color=gg_SB.3.3.mm),size=0.5) + 
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_2 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.2))+
  geom_point(aes(color=gg_MB.2),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 2") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_3 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.3))+
  geom_point(aes(color=gg_MB.3),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 3") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_4 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.4))+
  geom_point(aes(color=gg_MB.4),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 4") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_6 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.6))+
  geom_point(aes(color=gg_MB.6),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 6") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_8 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.8))+
  geom_point(aes(color=gg_MB.8),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 8") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_9 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.9))+
  geom_point(aes(color=gg_MB.9),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 9") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_12 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.12))+
  geom_point(aes(color=gg_MB.12),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 12") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")


sp_tf33 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_SB.3.3.mm))+
  geom_point(aes(color=gg_SB.3.3.mm),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf2 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.2))+
  geom_point(aes(color=gg_MB.2),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 2") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf3 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.3))+
  geom_point(aes(color=gg_MB.3),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 3") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf4 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.4))+
  geom_point(aes(color=gg_MB.4),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 4") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf6 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.6))+
  geom_point(aes(color=gg_MB.6),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 6") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf8 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.8))+
  geom_point(aes(color=gg_MB.8),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 8") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf9 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.9))+
  geom_point(aes(color=gg_MB.9),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 9") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf12 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.12))+
  geom_point(aes(color=gg_MB.12),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = xy.ratio) +
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 12") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

png(file='./OptimalSMS_rsfMRI/Documents/Figures/All16_FixedCoords_AttenuatedCorrelations_9p_3-4ratio.png',width=16,height=12,units='in',res=500)
ggarrange(sp_33, sp_tf33, sp_2, sp_tf2, sp_3, sp_tf3, sp_4, sp_tf4, sp_6, sp_tf6, sp_8, sp_tf8, sp_9, sp_tf9, sp_12, sp_tf12,
           ncol = 2, nrow = 8 )
dev.off()


######Plot with ratio of axes freely determined
sp_33 = ggplot(netdf_9p, aes(z_SB.2.mm,z_SB.3.3.mm))+
  geom_point(aes(color=gg_SB.3.3.mm),size=0.5) + 
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_2 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.2))+
  geom_point(aes(color=gg_MB.2),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 2") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_3 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.3))+
  geom_point(aes(color=gg_MB.3),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 3") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_4 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.4))+
  geom_point(aes(color=gg_MB.4),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 4") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_6 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.6))+
  geom_point(aes(color=gg_MB.6),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 6") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_8 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.8))+
  geom_point(aes(color=gg_MB.8),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 8") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_9 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.9))+
  geom_point(aes(color=gg_MB.9),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 9") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_12 = ggplot(netdf_9p, aes(z_SB.2.mm,z_MB.12))+
  geom_point(aes(color=gg_MB.12),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 12") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")


sp_tf33 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_SB.3.3.mm))+
  geom_point(aes(color=gg_SB.3.3.mm),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf2 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.2))+
  geom_point(aes(color=gg_MB.2),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 2") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf3 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.3))+
  geom_point(aes(color=gg_MB.3),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 3") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf4 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.4))+
  geom_point(aes(color=gg_MB.4),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 4") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf6 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.6))+
  geom_point(aes(color=gg_MB.6),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 6") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf8 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.8))+
  geom_point(aes(color=gg_MB.8),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 8") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf9 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.9))+
  geom_point(aes(color=gg_MB.9),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 9") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

sp_tf12 = ggplot(netdf_9p_tf, aes(z_SB.2.mm,z_MB.12))+
  geom_point(aes(color=gg_MB.12),size=0.5)+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "MB 12") +
  scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")

png(file='./OptimalSMS_rsfMRI/Documents/Figures/All16_FixedCoords_AttenuatedCorrelations_9p_FreeRatio.png',width=14,height=12,units='in',res=500)
ggarrange(sp_33, sp_tf33, sp_2, sp_tf2, sp_3, sp_tf3, sp_4, sp_tf4, sp_6, sp_tf6, sp_8, sp_tf8, sp_9, sp_tf9, sp_12, sp_tf12,
          ncol = 2, nrow = 8 )
dev.off()

# 
# # Next, create plots of Cohen's D
# # these are not currently used:
# mysample=sample(1:nrow(netdf_9p),10000)
# netdf_9p_sample = netdf_9p[mysample,]
# 
# sp_1 = ggplot(netdf_9p, aes(cd_SB.2.mm,cd_SB.3.3.mm))+geom_point(aes(color=gg_SB.3.3.mm),size=0.1)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = 1) +geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")+coord_fixed()
# 
# sp_2 = ggplot(netdf_9p, aes(cd_SB.2.mm,cd_MB.2))+geom_point(aes(color=gg_MB.2),size=0.1)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = 1) +geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")+coord_fixed()
# 
# sp_4 = ggplot(netdf_9p, aes(cd_SB.2.mm,cd_MB.4))+geom_point(aes(color=gg_MB.4),size=0.1)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = 1) +geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")+coord_fixed()
# 
# sp_8 = ggplot(netdf_9p, aes(cd_SB.2.mm,cd_MB.8))+geom_point(aes(color=gg_MB.8),size=0.1)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = 1) +geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")+coord_fixed()
# 
# sp_12 = ggplot(netdf_9p, aes(cd_SB.2.mm,cd_MB.12))+geom_point(aes(color=gg_MB.12),size=0.1)+geom_hline(yintercept=0)+geom_vline(xintercept=0)+coord_fixed(ratio = 1) +scale_color_gradientn(colours=pal,limits=c(0.9,4),oob = scales::squish,name = "g*")+geom_abline(slope=1,color='black')+ labs(x = "SB 2 mm", y = "SB 3.3 mm") +coord_fixed()
# 
# ggarrange(sp_1,sp_2,sp_4,sp_8,sp_12)
# # these don't seem very informative.
# 
# 
# 
# # tabulate the number of points in each quadrant:
# # not currently used:
# with(netdf_9p[netdf_9p$cd_MB.2>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.2>cd_SB.2.mm))
# with(netdf_9p[netdf_9p$cd_MB.3>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.3>cd_SB.2.mm))
# with(netdf_9p[netdf_9p$cd_MB.4>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.4>cd_SB.2.mm))
# with(netdf_9p[netdf_9p$cd_MB.8>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.8>cd_SB.2.mm))
# with(netdf_9p[netdf_9p$cd_MB.9>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.9>cd_SB.2.mm))
# with(netdf_9p[netdf_9p$cd_MB.12>0 & netdf_9p$cd_SB.2.mm>0,],mean(cd_MB.12>cd_SB.2.mm))
# 
# 
# 
# 
###############################################################
# Tabulate the number of significant correlations:
# note: package matlab has a function called sum!!! unload this function
#
###############################################################
detach(package:matlab, unload=TRUE)

#########


# Calculate number of significant edges using the number of samples (power may differ due to differences in number of scans)
  # create array of counts

  # old code if using the t-statistics. 
  nSubjectsArray = apply(allmb.cor.srt,c(1,2,3),function(x) sum(!is.na(x)))
  
  # convert Cohen's D to t-statistics
  tstat.z.cor = sqrt(nSubjectsArray)*cd.z.cor
  dim(tstat.z.cor)
  
  # check a few to ensure that are equal to t test:
    tstat.z.cor[3,1,1]
    t.test(allmb.z.cor.srt[3,1,1,])
    
    tstat.z.cor[230,235,8]
    t.test(allmb.z.cor.srt[230,235,8,])
    # looks good
  


###### Standardize number of subjects so power is equal:
# convert Cohen's D to t-statistics
  
# Exclude nodes that had some voxels with 0s in SB 3.3
power_resort = read.csv('./OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv',header = FALSE)
names(power_resort) = c('Community','ROI','BRisk_ROI')
exclude = power_resort[power_resort$ROI%in%c(4,5,9,10,250,83),3]
somemb.z.cor.srt = allmb.z.cor.srt[-exclude,-exclude,,]
somemb_tf.z.cor.srt = allmb_tf.z.cor.srt[-exclude,-exclude,,]
someLabelsLong = labelsLong[-exclude]

prop_sig_notf = calc.prop.active(somemb.z.cor.srt,p.thresh=0.05/choose(length(someLabelsLong),2),indexvariable=someLabelsLong,colLabels = varnames)
prop_sig_notf_long = pivot_longer(prop_sig_notf,cols=SB.3.3.mm:MB.12,names_to='Acquisition',values_to='Proportion')
prop_sig_notf_long$Acquisition = factor(prop_sig_notf_long$Acquisition,levels= c('SB.3.3.mm','SB.2.mm','MB.2','MB.3','MB.4','MB.6','MB.8','MB.9','MB.12'))
#prop_sig_notf_long$Community = factor(prop_sig_notf_long$Community,levels= prop_sig_notf$Community)

prop_sig_tf = calc.prop.active(somemb_tf.z.cor.srt,p.thresh=0.05/choose(length(someLabelsLong),2),indexvariable=someLabelsLong,colLabels = varnames)
prop_sig_tf_long = pivot_longer(prop_sig_tf,cols=SB.3.3.mm:MB.12,names_to='Acquisition',values_to='Proportion')
prop_sig_tf_long$Acquisition = factor(prop_sig_tf_long$Acquisition,levels= c('SB.3.3.mm','SB.2.mm','MB.2','MB.3','MB.4','MB.6','MB.8','MB.9','MB.12'))

gg1 = ggplot(data = prop_sig_notf_long[prop_sig_notf_long$Community!="Uncertain",], aes(x=Community, y=Proportion, fill=Acquisition))+geom_bar(stat='identity',pos='dodge',colour='darkslategrey')+ylab('Proportion Edges Activated')+ylim(c(0,0.18))+ggtitle("A) 9p")+theme(axis.text=element_text(size=12))+scale_fill_brewer(palette="Paired")

gg2 = ggplot(data = prop_sig_tf_long[prop_sig_tf_long$Community!="Uncertain",], aes(x=Community, y=Proportion, fill=Acquisition))+geom_bar(stat='identity',pos='dodge',colour='darkslategrey')+ylab('Proportion Edges Activated')+ylim(c(0,0.18))+ggtitle("B) 9p+bandpass")+theme(axis.text=element_text(size=12))+scale_fill_brewer(palette = "Paired")

png(file='./OptimalSMS_rsfMRI/Documents/Figures/ProportionActivated_bonf.png',width=14,height=12,units='in',res=500)
ggarrange(gg1, gg2, nrow=2 )
dev.off()

### Calculate mean cohen's d across all edges in a community
cd_notf=calc.cohensd.comm(somemb.z.cor.srt,indexvariable = someLabelsLong,colLabels = varnames)
cd_notf_long=pivot_longer(cd_notf,cols=SB.3.3.mm:MB.12,names_to='Acquisition',values_to='CohensD')
cd_notf_long$Acquisition = factor(cd_notf_long$Acquisition,levels= c('SB.3.3.mm','SB.2.mm','MB.2','MB.3','MB.4','MB.6','MB.8','MB.9','MB.12'))

cd_tf=calc.cohensd.comm(somemb_tf.z.cor.srt,indexvariable = someLabelsLong,colLabels = varnames)
cd_tf_long=pivot_longer(cd_tf,cols=SB.3.3.mm:MB.12,names_to='Acquisition',values_to='CohensD')
cd_tf_long$Acquisition = factor(cd_tf_long$Acquisition,levels= c('SB.3.3.mm','SB.2.mm','MB.2','MB.3','MB.4','MB.6','MB.8','MB.9','MB.12'))

gg1 = ggplot(data = cd_notf_long[cd_notf_long$Community!="Uncertain",], aes(x=Community, y=CohensD, fill=Acquisition))+geom_bar(stat='identity',pos='dodge',colour='darkslategrey')+ylab('Mean Cohen\'s d')+ggtitle("A) 9p")+theme(axis.text=element_text(size=12))+ylim(c(0,0.7))+scale_fill_brewer(palette="Paired")

gg2 = ggplot(data = cd_tf_long[cd_tf_long$Community!="Uncertain",], aes(x=Community, y=CohensD, fill=Acquisition))+geom_bar(stat='identity',pos='dodge',colour='darkslategrey')+ylab('Mean Cohen\'s d')+ylim(c(0,0.7))+ggtitle("B) 9p+bandpass")+theme(axis.text=element_text(size=12))+scale_fill_brewer(palette="Paired")

png(file='./OptimalSMS_rsfMRI/Documents/Figures/MeanCohensD.png',width=14,height=12,units='in',res=500)
ggarrange(gg1, gg2, nrow=2 )
dev.off()


save(somemb.z.cor.srt,somemb_tf.z.cor.srt,someLabelsLong,varnames,prop_sig_tf,prop_sig_notf,cd_notf,cd_tf,file='~/risk_share/OptimalSMS_rsfMRI/Results/somemb.z.cor.srt.RData')
