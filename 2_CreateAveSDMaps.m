%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin Risk
% Calculate average SD Maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('~/mfunctions/NIfTI_20140122')
% resultspath = '~/risk_share/OptimalSMS_rsfMRI/Results/';
% datapath = '~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/';

mblist={'rsfMRI_3_3mm_iso_MB1_AP','rsfMRI_2mm_iso_MB1_AP','rsfMRI_2mm_iso_MB2_AP','rsfMRI_2mm_iso_MB3_AP','rsfMRI_2mm_iso_MB4_AP',...
    'rsfMRI_2mm_iso_MB6_AP','rsfMRI_2mm_iso_MB8_AP', 'rsfMRI_2mm_iso_MB9_AP', 'rsfMRI_2mm_iso_MB12_AP'};

% Note: subjects must be in the same order as gross motion Subjects
excludeScans = readmatrix('~/risk_share/OptimalSMS_rsfMRI/Results/exclude_subjects.csv');

subjectList=excludeScans(:,1);

nSubject = length(subjectList);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create nifti files with average standard deviations across subjects for
% each MB:
stdev_9p_Array = NaN(91,109,91,9,32);
stdev_9p_tf_Array = NaN(91,109,91,9,32);
% stdev_aroma_Array = NaN(91,109,91,9,32);
stdev_9p_sm_Array = NaN(91,109,91,9,32);
stdev_9p_tf_sm_Array = NaN(91,109,91,9,32);
% grab header etc from an example subject

for j=1:length(subjectList)
    subject = num2str(subjectList(j));
    for i = 1:length(mblist)
      flag = excludeScans(j,i+1);
      if flag==0
          try
              im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subject '/nii/' mblist{i} '/' mblist{i} '_pp_9p_stdev.nii.gz'];
              stdev_9p_Array(:,:,:,i,j)=niftiread(im_file);

              im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subject '/nii/' mblist{i} '/' mblist{i} '_pp_9p_tf_stdev.nii.gz'];
              stdev_9p_tf_Array(:,:,:,i,j)=niftiread(im_file);

  %            im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subject '/nii/' mblist{i} '/' mblist{i} '_pp_9p_aroma_stdev.nii.gz'];
  %            stdev_aroma_Array(:,:,:,i,j)=niftiread(im_file);
              
              im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subject '/nii/' mblist{i} '/' mblist{i} '_pp_9p_sm_stdev.nii.gz'];
              stdev_9p_sm_Array(:,:,:,i,j)=niftiread(im_file);

              im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subject '/nii/' mblist{i} '/' mblist{i} '_pp_9p_tf_sm_stdev.nii.gz'];
              stdev_9p_tf_sm_Array(:,:,:,i,j)=niftiread(im_file);


          catch
               warning(['No ' mblist{i} ' for participant ' subject])
          end
      end
    end
end

mean_stdev_9p=mean(stdev_9p_Array,5,'omitnan');
mean_stdev_9p_tf=mean(stdev_9p_tf_Array,5,'omitnan');
% mean_stdev_aroma=mean(stdev_aroma_Array,5,'omitnan');
mean_stdev_9p_sm=mean(stdev_9p_sm_Array,5,'omitnan');
mean_stdev_9p_tf_sm=mean(stdev_9p_tf_sm_Array,5,'omitnan');

tempinfo = niftiinfo('~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/102954/nii/rsfMRI_2mm_iso_MB1_AP/rsfMRI_2mm_iso_MB1_AP_pp_9p_tf.nii.gz');
tempinfo.ImageSize=[91 109 91 9];

niftiwrite(single(mean_stdev_9p),'~/risk_share/OptimalSMS_rsfMRI/Results/mean_stdev_9p.nii',tempinfo);
niftiwrite(single(mean_stdev_9p_tf),'~/risk_share/OptimalSMS_rsfMRI/Results/mean_stdev_9p_tf.nii',tempinfo);
% niftiwrite(single(mean_stdev_aroma),'~/risk_share/OptimalSMS_rsfMRI/Results/mean_stdev_9p_aroma.nii',tempinfo);
niftiwrite(single(mean_stdev_9p_sm),'~/risk_share/OptimalSMS_rsfMRI/Results/mean_stdev_9p_sm.nii',tempinfo);
niftiwrite(single(mean_stdev_9p_tf_sm),'~/risk_share/OptimalSMS_rsfMRI/Results/mean_stdev_9p_tf_sm.nii',tempinfo);




