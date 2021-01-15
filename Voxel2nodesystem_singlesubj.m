% BRisk
%
% Calculate average time series and correlation matrices for multiple preprocessing pipelines:
%   1. 9P with no temporal filtering
%   2. 9P with temporal filtering
%   3. 9p with no temporal filtering and with 6-mm smoothing
%   4. 9p with temporal filtering and with 6-mm smoothing
%   5. ICA AROMA with no temporal filtering

addpath('~/mfunctions/NIfTI_20140122')
subj='subjectID';
s=load_nii('~/risk_share/OptimalSMS_rsfMRI/PowerAtlas/total_sphere.nii');  

nNode=264;
mblist={'rsfMRI_3_3mm_iso_MB1_AP','rsfMRI_2mm_iso_MB1_AP','rsfMRI_2mm_iso_MB2_AP','rsfMRI_2mm_iso_MB3_AP','rsfMRI_2mm_iso_MB4_AP',...
    'rsfMRI_2mm_iso_MB6_AP','rsfMRI_2mm_iso_MB8_AP', 'rsfMRI_2mm_iso_MB9_AP', 'rsfMRI_2mm_iso_MB12_AP'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. 9P Preprocessing with no temporal filtering:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parfor j=1:length(subjList) 
% parfor does not work on outside of multiple loops, must be furthest in

for i = 1:length(mblist)
  try
  im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{i} '/' mblist{i} '_pp_9p.nii.gz'];
  s2=load_nii(im_file);
  [nx,ny,nz,nt] = size(s2.img);
  node_ts = zeros(nt,nNode);
    for n=1:nNode %loop through nodes
            n_ind=find(s.img==n);
            myimg = reshape(s2.img,[nx*ny*nz,nt]);
            sub_myimg = myimg(n_ind,:);
            node_ts(:,n)=mean(sub_myimg,1);
    end

    ot_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/' mblist{i} '_pp_9p_pwrts.mat'];
    corrmat = corr(node_ts);
    sdvec = std(node_ts)';
    save(ot_file,'node_ts','corrmat','sdvec');
    fprintf(['Finished ' mblist{i} '\n*********\n'])
  catch
            warning(['No ' mblist{i} '_9p for participant ' subj])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. 9P Preprocessing with temporal filtering (9p+bp):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(mblist)
  try
  im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{i} '/' mblist{i} '_pp_9p_tf.nii.gz'];
  s2=load_nii(im_file);
  [nx,ny,nz,nt] = size(s2.img);
  node_ts = zeros(nt,nNode);
    for n=1:nNode %loop through nodes
            n_ind=find(s.img==n);
            myimg = reshape(s2.img,[nx*ny*nz,nt]);
            sub_myimg = myimg(n_ind,:);
            node_ts(:,n)=mean(sub_myimg,1);
    end

    ot_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/' mblist{i} '_pp_9p_tf_pwrts.mat'];
    corrmat = corr(node_ts);
    sdvec = std(node_ts)';
    save(ot_file,'node_ts','corrmat','sdvec');
    fprintf(['Finished ' mblist{i} '\n*********\n'])
  catch
            warning(['No ' mblist{i} '_9p_tf for participant ' subj])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. 9P Preprocessing with spatial smoothing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(mblist)
  try
  im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{i} '/' mblist{i} '_pp_9p_sm.nii.gz'];
  s2=load_nii(im_file);
  [nx,ny,nz,nt] = size(s2.img);
  node_ts = zeros(nt,nNode);
    for n=1:nNode %loop through nodes
            n_ind=find(s.img==n);
            myimg = reshape(s2.img,[nx*ny*nz,nt]);
            sub_myimg = myimg(n_ind,:);
            node_ts(:,n)=mean(sub_myimg,1);
    end

    ot_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/' mblist{i} '_pp_9p_sm_pwrts.mat'];
    corrmat = corr(node_ts);
    sdvec = std(node_ts)';
    save(ot_file,'node_ts','corrmat','sdvec');
    fprintf(['Finished ' mblist{i} '\n*********\n'])
  catch
            warning(['No ' mblist{i} '9p_sm for participant ' subj])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. 9P Preprocessing with temporal filtering and spatial smoothing:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(mblist)
  try
  im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{i} '/' mblist{i} '_pp_9p_tf_sm.nii.gz'];
  s2=load_nii(im_file);
  [nx,ny,nz,nt] = size(s2.img);
  node_ts = zeros(nt,nNode);
    for n=1:nNode %loop through nodes
            n_ind=find(s.img==n);
            myimg = reshape(s2.img,[nx*ny*nz,nt]);
            sub_myimg = myimg(n_ind,:);
            node_ts(:,n)=mean(sub_myimg,1);
    end

    ot_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/' mblist{i} '_pp_9p_tf_sm_pwrts.mat'];
    corrmat = corr(node_ts);
    sdvec = std(node_ts)';
    save(ot_file,'node_ts','corrmat','sdvec');
    fprintf(['Finished ' mblist{i} '\n*********\n'])
  catch
            warning(['No ' mblist{i} '9p_tf_sm for participant ' subj])
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. 9P Preprocessing with ICA AROMA:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for i = 1:length(mblist)
%   try
%   im_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{i} '/' mblist{i} '_pp_9p_aroma.nii.gz'];
%   s2=load_nii(im_file);
%   [nx,ny,nz,nt] = size(s2.img);
%   node_ts = zeros(nt,nNode);
%     for n=1:nNode %loop through nodes
%             n_ind=find(s.img==n);
%             myimg = reshape(s2.img,[nx*ny*nz,nt]);
%             sub_myimg = myimg(n_ind,:);
%             node_ts(:,n)=mean(sub_myimg,1);
%     end
% 
%     ot_file=['~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/' mblist{i} '_pp_9p_aroma_pwrts.mat'];
%     corrmat = corr(node_ts);
%     sdvec = std(node_ts)';
%     save(ot_file,'node_ts','corrmat','sdvec');
%     fprintf(['Finished ' mblist{i} '\n*********\n'])
%   catch
%             warning(['No ' mblist{i} '9p_aroma for participant ' subj])
%   end
% end

