########################
# Benjamin Risk
# Examine motion
# Grab motion files and create a 
# concatenated data frame

# get subject list:
subjectList = dir('~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/')

# remove the extra files:

subjectList = subjectList[1:32]

subjectList

mbList=c("rsfMRI_3_3mm_iso_MB1_AP","rsfMRI_2mm_iso_MB1_AP","rsfMRI_2mm_iso_MB2_AP","rsfMRI_2mm_iso_MB3_AP","rsfMRI_2mm_iso_MB4_AP","rsfMRI_2mm_iso_MB6_AP","rsfMRI_2mm_iso_MB8_AP","rsfMRI_2mm_iso_MB9_AP","rsfMRI_2mm_iso_MB12_AP")


meanrms_motion=data.frame('Subject'=subjectList,'SB3p3'=0,'SB2'=0,'MB2'=0,'MB3'=0,'MB4'=0,'MB6'=0,'MB8'=0,'MB9'=0,'MB12'=0)

ngr25_motion = meanrms_motion

# Ciric criteria: 
#   a. mean relative RMS >0.2 mm
#   b. 16% exceeded 0.25
for (i in 1:32) {
  subj=subjectList[i]
  k=1
  setwd(paste0('~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/',subj,'/nii/'))
  for (j in mbList) {
    k=k+1
    temp = try(read.table(paste0('./',j,'/',j,'_mcf.nii.gz_rel.rms')))
    if(class(temp)!='try-error') {
      meanrms_motion[i,k] = mean(temp[,1])
      ngr25_motion[i,k] = mean(temp[,1]>0.25)
    }
  }
}

# 4 volumes with meanrms that fail criteria:
sum(meanrms_motion[,2:10]>0.2)

# 9 fail due to >0.25 mm fd in 1/6 of volumes:
sum(ngr25_motion[,2:10]>(1/6))


exclude = (meanrms_motion[,2:10]>0.2)|(ngr25_motion[,2:10]>(1/6))

sum(exclude)
# 9 scans will be removed


exclude_scans = data.frame("Subjects"=subjectList,exclude)

# exclude subject 193728 MB9, which had the wrong TR
# also "exclude" 195590 MB9, which was not collected`
exclude_scans[subjectList%in%c('193728','195590'),c('MB9')]=TRUE

exclude_scans

apply(exclude_scans[,c(2:10)],2,table)

save(file='~/risk_share/OptimalSMS_rsfMRI/Results/GrossMotion.RData',meanrms_motion,ngr25_motion,exclude_scans)

exclude_scans[,2:10] = 1*(exclude_scans[,2:10])

write.csv(exclude_scans,file = '~/risk_share/OptimalSMS_rsfMRI/Results/exclude_subjects.csv',row.names=FALSE)
write.csv(meanrms_motion,file = '~/risk_share/OptimalSMS_rsfMRI/Results/meanrms_motion.csv',row.names=FALSE)
write.csv(ngr25_motion,file='~/risk_share/OptimalSMS_rsfMRI/Results/ngr25_motion.csv',row.names=FALSE)

######
# Subject with median motion for use in plots

meandf = data.frame('Subjects'=exclude_scans$Subjects,'MeanRMS'=apply(meanrms_motion[,2:10],1,mean))
meandf$meanmean=mean(meandf[,2])
meandf$diff=abs(meandf$MeanRMS-meandf$meanmean)
meandf
meandf$Subject[which.min(meandf$diff)]




