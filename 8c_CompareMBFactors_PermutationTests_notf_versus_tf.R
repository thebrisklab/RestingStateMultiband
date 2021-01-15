
#################################
# Benjamin Risk
# September 8 2020
# Permutation test for difference in activation
# with and without temporal filtering
####################################
a = .libPaths()
.libPaths(c('/home/benjamin.risk/Rlibs',a))

.libPaths()


library(doParallel)
setwd('~/risk_share')
source('./OptimalSMS_rsfMRI/Programs/Functions/fconn_fun.R')
load('./OptimalSMS_rsfMRI/Results/somemb.z.cor.srt.RData')
subj_exclude=read.csv('./OptimalSMS_rsfMRI/Results/exclude_subjects.csv')

# note: this dataset already has subjects excluded
# subj_excluded used when subsetting to subjects that match between datasets

# acquisition indices:
# 1: SB 3.3
# 2: SB 2 mm
# 3: MB 2
# 4: MB 3
# 5: MB 4
# 6: MB 6
# 7: MB 8
# 8: MB 9
# 9: MB 12

no_cores = detectCores()-11
registerDoParallel(no_cores)
seed.vector = c(1:10000)
nnode=dim(somemb_tf.z.cor.srt)[1]
p.thresh = 0.05/choose(nnode,2)


a = proc.time()
sb3.3 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=1,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
sb2 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=2,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb2 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=3,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb3 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=4,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb4 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=5,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb6 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=6,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb8 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb9 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=8,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb12 = perm.test.prop.activated.notf_v_tf(datain_notf=somemb.z.cor.srt,datain_tf=somemb_tf.z.cor.srt,subject_exclude_df=subj_exclude,acquisition=9,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)

time = proc.time()-a

p.array.notf.tf = matrix(NA,15,9)
colnames(p.array.notf.tf) = varnames
rownames(p.array.notf.tf) = names(sb2[[2]])

p.array.notf.tf[,1] = sb3.3$p.values
p.array.notf.tf[,2] = sb2$p.values
p.array.notf.tf[,3] = mb2$p.values
p.array.notf.tf[,4] = mb3$p.values
p.array.notf.tf[,5] = mb4$p.values
p.array.notf.tf[,6] = mb6$p.values
p.array.notf.tf[,7] = mb8$p.values
p.array.notf.tf[,8] = mb9$p.values
p.array.notf.tf[,9] = mb12$p.values

save(file='~/risk_share/OptimalSMS_rsfMRI/Results/permtests_edge_density_notf_versus_tf.RData',p.array.notf.tf,time,sb3.3,sb2,mb2,mb3,mb4,mb6,mb8,mb9,mb12)
