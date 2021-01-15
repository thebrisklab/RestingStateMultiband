addpath('~/risk_jhu_share/OptimalSMS_rsfMRI/mfunctions/NIfTI_20140122')

SeedMaps_create({'_pp_9p', '_pp_9p_tf','_pp_9p_sm','_pp_9p_tf_sm'}, 'DorsalRostralPutamenSphere')
SeedMaps_create({'_pp_9p', '_pp_9p_tf','_pp_9p_sm','_pp_9p_tf_sm'}, 'LM1sphere')
% takes about 20 minutes per pipeline per roi


% There are five sets of files:
%
%rsfMRI_2mm_iso_MB1_AP_pp_9p.nii.gz The data processed without bandpass
%filtering

%rsfMRI_2mm_iso_MB1_AP_pp_9p_tf.nii.gz The data processed with bandpass
%filtering

%rsfMRI_2mm_iso_MB1_AP_pp_9p_sm.nii.gz The data processed without bandpass
%filtering and with 6-mm FWHM smoothing

%rsfMRI_2mm_iso_MB1_AP_pp_9p_tf.nii.gz The data processed with bandpass
%filtering and 6-mm FWHM smoothing

% There is also rsfMRI_2mm_iso_MB1_AP_pp_aroma.nii.gz. Initial
% inspection indicated AROMA remove a large amount of signal,
% and detailed analyses are not being conducted with this pipeline.
