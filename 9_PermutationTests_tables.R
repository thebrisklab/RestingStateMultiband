#######################
# Benjamin Risk
# Create tables from permutation tests of differences between edge density for MB 8 and other acquisitions
# Create table of differences between no bandpass and bandpass
# Modify this depending on whether running as an interactive job on the 
# cluster versus locally:

setwd('~/Dropbox/OptimalSMS_rsfMRI/Results')

#load('somemb.z.cor.srt.RData')
library(xtable)
load('permtests_edge_density_notf.RData')

diff.for.xtable = function(Rdata,alpha.level=0.01) {
  load(Rdata)
  diff.array = matrix(NA,15,8)
  #colnames(diff.array) = varnames[c(1:6,8:9)]
  rownames(diff.array) = names(sb2.mb8[[2]])
  diff.array[,1] = sb3.3.mb8$observed.diff
  diff.array[,2] = sb2.mb8$observed.diff
  diff.array[,3] = mb2.mb8$observed.diff
  diff.array[,4] = mb3.mb8$observed.diff
  diff.array[,5] = mb4.mb8$observed.diff
  diff.array[,6] = mb6.mb8$observed.diff
  diff.array[,7] = mb9.mb8$observed.diff
  diff.array[,8] = mb12.mb8$observed.diff
  
  diff.pvalue=matrix(NA,15,8)
  for (j in 1:8) {
    for (i in 1:15) {
      if (p.array[i,j]<alpha.level) {
        diff.pvalue[i,j]=paste0('\\bf{',round(diff.array[i,j],3)," p=",round(p.array[i,j],3),'}')
      }  
      else {
      diff.pvalue[i,j]=paste0(round(diff.array[i,j],3)," p=",round(p.array[i,j],3))
      }
    }
  }
  colnames(diff.pvalue) = c('SB 3.3 mm','SB 2 mm','MB 2', 'MB 3', 'MB 4', 'MB 6','MB 9','MB 12')
  rownames(diff.pvalue) = rownames(p.array)
  diff.pvalue = diff.pvalue[-12,]

  # use longer version of label names:
  labels=c('Auditory','Cerebellum','CO-Task Ctrl','Default Mode','Dorsal Att','FP-Task Ctrl','Memory','Salience','Som-Hand','Som-Mouth','Subcortical','Ventral Att','Visual','All')
  rownames(diff.pvalue)=labels
  diff.pvalue
}

diff.pvalue.notf = diff.for.xtable('permtests_edge_density_notf.RData')

print(xtable(diff.pvalue.notf[,1:4]),type='latex',sanitize.text.function = function(x) {x})
print(xtable(diff.pvalue.notf[,5:8]),type='latex',sanitize.text.function = function(x) {x})

## Create tables for 9p+bandpass:

load('permtests_edge_density_tf.RData')
diff.pvalue.tf = diff.for.xtable('permtests_edge_density_tf.RData')
print(xtable(diff.pvalue.tf[,1:4]),type='latex',sanitize.text.function = function(x) {x})
print(xtable(diff.pvalue.tf[,5:8]),type='latex',sanitize.text.function = function(x) {x})


#########################


load('permtests_edge_density_notf_versus_tf.RData')
diff.array.notf_vs_tf = matrix(NA,15,9)
rownames(diff.array.notf_vs_tf) = names(sb2[[2]])
diff.array.notf_vs_tf[,1] = sb3.3$observed.diff
diff.array.notf_vs_tf[,2] = sb2$observed.diff
diff.array.notf_vs_tf[,3] = mb2$observed.diff
diff.array.notf_vs_tf[,4] = mb3$observed.diff
diff.array.notf_vs_tf[,5] = mb4$observed.diff
diff.array.notf_vs_tf[,6] = mb6$observed.diff
diff.array.notf_vs_tf[,7] = mb8$observed.diff
diff.array.notf_vs_tf[,8] = mb9$observed.diff
diff.array.notf_vs_tf[,9] = mb12$observed.diff


diff.pvalue.notf_vs_tf=matrix(NA,15,9)
alpha.level=0.01
for (j in 1:9) {
  for (i in 1:15) {
    if (p.array.notf.tf[i,j]<alpha.level) {
      diff.pvalue.notf_vs_tf[i,j]=paste0('\\bf{',round(diff.array.notf_vs_tf[i,j],3)," p=",round(p.array.notf.tf[i,j],3),'}')
    }  
    else {
      diff.pvalue.notf_vs_tf[i,j]=paste0(round(diff.array.notf_vs_tf[i,j],3)," p=",round(p.array.notf.tf[i,j],3))
    }
  }
}

colnames(diff.pvalue.notf_vs_tf)=varnames
rownames(diff.pvalue.notf_vs_tf)=names(mb12$observed.diff)

diff.pvalue.notf_vs_tf
diff.pvalue.notf_vs_tf = diff.pvalue.notf_vs_tf[-12,]
rownames(diff.pvalue.notf_vs_tf)= c('Auditory','Cerebellum','CO-Task Ctrl','Default Mode','Dorsal Att','FP-Task Ctrl','Memory','Salience','Som-Hand','Som-Mouth','Subcortical','Ventral Att','Visual','All')

# start here
print(xtable(diff.pvalue.notf_vs_tf[,1:5]),type='latex',sanitize.text.function=function(x) {x})
print(xtable(diff.pvalue.notf_vs_tf[,6:9]),type='latex',sanitize.text.function=function(x) {x})

######################################
# Create tables for smoothed data:
load('permtests_edge_density_notf.RData')
diff.pvalue.notf.sm = diff.for.xtable('permtests_edge_density_notf_sm.RData')

print(xtable(diff.pvalue.notf.sm[,1:4]),type='latex',sanitize.text.function = function(x) {x})
print(xtable(diff.pvalue.notf.sm[,5:8]),type='latex',sanitize.text.function = function(x) {x})


## NOT INCLUDED IN MANUSCRIPT OR SUPPLEMENT:
load('permtests_edge_density_tf.RData')
diff.pvalue.tf.sm = diff.for.xtable('permtests_edge_density_tf_sm.RData')
print(xtable(diff.pvalue.tf.sm[,1:4]),type='latex',sanitize.text.function = function(x) {x})
print(xtable(diff.pvalue.tf.sm[,5:8]),type='latex',sanitize.text.function = function(x) {x})
