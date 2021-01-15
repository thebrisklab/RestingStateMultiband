
#################################
# Benjamin Risk
# August 23 2020
# Permutation test for difference in activation:
####################################
a = .libPaths()
.libPaths(c('/home/benjamin.risk/Rlibs',a))

.libPaths()


library(doParallel)
setwd('~/risk_share')
source('./OptimalSMS_rsfMRI/Programs/Functions/fconn_fun.R')
load('./OptimalSMS_rsfMRI/Results/somemb.z.cor.srt.RData')
subj_exclude=read.csv('~/risk_share/OptimalSMS_rsfMRI/Results/exclude_subjects.csv')

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
nnode=dim(somemb.z.cor.srt)[1]
p.thresh = 0.05/choose(nnode,2)

a = proc.time()
sb3.3.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=1,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
sb2.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=2,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb2.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=3,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb3.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=4,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb4.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=5,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb6.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=6,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb9.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=8,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
mb12.mb8 = perm.test.prop.activated(datain=somemb.z.cor.srt,subject_exclude_df=subj_exclude,subject_match=TRUE,acquisition1=9,acquisition2=7,p.threshold=p.thresh,indexvariable=someLabelsLong,seed.vector=seed.vector)
time = proc.time()-a

p.array = matrix(NA,15,8)
colnames(p.array) = varnames[c(1:6,8:9)]
rownames(p.array) = names(sb2.mb8[[2]])

p.array[,1] = sb3.3.mb8$p.values
p.array[,2] = sb2.mb8$p.values
p.array[,3] = mb2.mb8$p.values
p.array[,4] = mb3.mb8$p.values
p.array[,5] = mb4.mb8$p.values
p.array[,6] = mb6.mb8$p.values
p.array[,7] = mb9.mb8$p.values
p.array[,8] = mb12.mb8$p.values

save(file='~/risk_share/OptimalSMS_rsfMRI/Results/permtests_edge_density_notf.RData',p.array,time,sb2.mb8,sb3.3.mb8,mb2.mb8,mb3.mb8,mb4.mb8,mb6.mb8,mb9.mb8,mb12.mb8)
