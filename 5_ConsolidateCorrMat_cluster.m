%%%%%%%%%%%%%%%%%%%%%%%%%%
% BRisk
% Calculate average time courses and the correlation matrices
%   for the Power atlas. 
%   NA subjects that fail the Ciric threshold
% Uses the output from Voxel2nodesystem_loopSubjects.sh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/mfunctions/NIfTI_20140122')
resultspath = '~/risk_share/OptimalSMS_rsfMRI/Results/';
datapath = '~/risk_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/';

mblist={'rsfMRI_3_3mm_iso_MB1_AP','rsfMRI_2mm_iso_MB1_AP','rsfMRI_2mm_iso_MB2_AP','rsfMRI_2mm_iso_MB3_AP','rsfMRI_2mm_iso_MB4_AP',...
    'rsfMRI_2mm_iso_MB6_AP','rsfMRI_2mm_iso_MB8_AP', 'rsfMRI_2mm_iso_MB9_AP', 'rsfMRI_2mm_iso_MB12_AP'};

excludeScans = readmatrix('~/risk_share/OptimalSMS_rsfMRI/Results/exclude_subjects.csv');

subjectList=excludeScans(:,1);

nSubject = length(subjectList);

nNode=264;
allmb_cor = NaN(nNode,nNode,9,32);
allmb_sdsd = NaN(nNode,nNode,9,32);
allmb_gg = NaN(nNode,nNode,9,32);
allmb_sd_node = NaN(nNode,9,32);

allmb_tf_cor = NaN(nNode,nNode,9,32);
allmb_tf_sdsd = NaN(nNode,nNode,9,32);
allmb_tf_gg = NaN(nNode,nNode,9,32);
allmb_tf_sd_node = NaN(nNode,9,32);

allmb_sm_cor = NaN(nNode,nNode,9,32);
allmb_sm_sdsd = NaN(nNode,nNode,9,32);
allmb_sm_gg = NaN(nNode,nNode,9,32);
allmb_sm_sd_node = NaN(nNode,9,32);

allmb_sm_tf_cor = NaN(nNode,nNode,9,32);
allmb_sm_tf_sdsd = NaN(nNode,nNode,9,32);
allmb_sm_tf_gg = NaN(nNode,nNode,9,32);
allmb_sm_tf_sd_node = NaN(nNode,9,32);

% allmb_aroma_cor = NaN(nNode,nNode,9,32);
% allmb_aroma_sdsd = NaN(nNode,nNode,9,32);
% allmb_aroma_gg = NaN(nNode,nNode,9,32);
% allmb_aroma_sd_node = NaN(nNode,9,32);

tic;
for j=1:length(subjectList)
    subject = num2str(subjectList(j));
    for i=1:9 
        flag = excludeScans(j,i+1);
        if flag==0
            try
                load([datapath subject '/' mblist{i} '_pp_9p_pwrts.mat']);
                allmb_cor(:,:,i,j) = corrmat;
                allmb_sdsd(:,:,i,j) = sdvec*sdvec';
                allmb_sd_node(:,i,j) = sdvec;
                
                load([datapath subject '/' mblist{i} '_pp_9p_tf_pwrts.mat']);
                allmb_tf_cor(:,:,i,j) = corrmat;
                allmb_tf_sdsd(:,:,i,j) = sdvec*sdvec';
                allmb_tf_sd_node(:,i,j) = sdvec;
                
                load([datapath subject '/' mblist{i} '_pp_9p_sm_pwrts.mat']);
                allmb_sm_cor(:,:,i,j) = corrmat;
                allmb_sm_sdsd(:,:,i,j) = sdvec*sdvec';
                allmb_sm_sd_node(:,i,j) = sdvec;
                
                load([datapath subject '/' mblist{i} '_pp_9p_tf_sm_pwrts.mat']);
                allmb_tf_sm_cor(:,:,i,j) = corrmat;
                allmb_tf_sm_sdsd(:,:,i,j) = sdvec*sdvec';
                allmb_tf_sm_sd_node(:,i,j) = sdvec;
                
          %      load([datapath subject '/' mblist{i} '_pp_9p_aroma_pwrts.mat']);
          %      allmb_aroma_cor(:,:,i,j) = corrmat;
          %      allmb_aroma_sdsd(:,:,i,j) = sdvec*sdvec';
          %      allmb_aroma_sd_node(:,i,j) = sdvec;
                
            catch
                warning(['File '  mblist{i} ' not found'])
            end
        end
    end

    for i=1:9
        allmb_gg(:,:,i,j) = allmb_sdsd(:,:,i,j)./allmb_sdsd(:,:,2,j);
        allmb_tf_gg(:,:,i,j) = allmb_tf_sdsd(:,:,i,j)./allmb_tf_sdsd(:,:,2,j);
        allmb_sm_gg(:,:,i,j) = allmb_sm_sdsd(:,:,i,j)./allmb_sm_sdsd(:,:,2,j);
        allmb_tf_sm_gg(:,:,i,j) = allmb_tf_sm_sdsd(:,:,i,j)./allmb_tf_sm_sdsd(:,:,2,j);
      %  allmb_aroma_gg(:,:,i,j) = allmb_aroma_sdsd(:,:,i,j)./allmb_aroma_sdsd(:,:,2,j);
    end
end
toc

% some subcortical areas are not completely covered by voxels in the
% MNI aligned space. Tabulate number of missings:

poweratlas = csvread('~/risk_share/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num.csv');
temp = squeeze(min(isnan(allmb_cor),[],2));
temp = squeeze(sum(temp,3));


fprintf('Create a table to look at whether all nodes contain data for all subjects:')
[poweratlas,temp]
fprintf('Note: In SB 3.3, node 10 is missing data for most subjects. Node 4 is also missing data. These \n are in the "uncertain" community.\n')


% Re-sort to put subcortical first:
% -1 is Uncertain -- change to 14
% subcortical is 10 -- change to -1
poweratlas2=poweratlas;
poweratlas2(poweratlas2(:,2)==-1,2)=14;
poweratlas2(poweratlas2(:,2)==10,2)=-1;
poweratlas2(poweratlas2(:,2)==13,2)=0;
tabulate(poweratlas2(:,2))

[sortedpower,indexpower] = sort(poweratlas2(:,2));

[sortedpower,indexpower,[1:264]']
% 240 corresponds to 4
% 246 correspond to 10
csvwrite('~/risk_share/OptimalSMS_rsfMRI/PowerAtlas/Neuron_consensus_264_abridged_num_resort.csv',[sortedpower,indexpower,[1:264]'])

% re-sort to put subcortical first:
allmb_cor_srt = allmb_cor(indexpower,indexpower,:,:);
allmb_sdsd_srt = allmb_sdsd(indexpower,indexpower,:,:);
allmb_sd_node_srt = allmb_sd_node(indexpower,:,:);
allmb_gg_srt = allmb_gg(indexpower,indexpower,:,:);


allmb_tf_cor_srt = allmb_tf_cor(indexpower,indexpower,:,:);
allmb_tf_sdsd_srt = allmb_tf_sdsd(indexpower,indexpower,:,:);
allmb_tf_gg_srt = allmb_tf_gg(indexpower,indexpower,:,:);
allmb_tf_sd_node_srt = allmb_tf_sd_node(indexpower,:,:);

allmb_sm_cor_srt = allmb_sm_cor(indexpower,indexpower,:,:);
allmb_sm_sdsd_srt = allmb_sm_sdsd(indexpower,indexpower,:,:);
allmb_sm_sd_node_srt = allmb_sm_sd_node(indexpower,:,:);
allmb_sm_gg_srt = allmb_sm_gg(indexpower,indexpower,:,:);


allmb_tf_sm_cor_srt = allmb_tf_sm_cor(indexpower,indexpower,:,:);
allmb_tf_sm_sdsd_srt = allmb_tf_sm_sdsd(indexpower,indexpower,:,:);
allmb_tf_sm_gg_srt = allmb_tf_sm_gg(indexpower,indexpower,:,:);
allmb_tf_sm_sd_node_srt = allmb_tf_sm_sd_node(indexpower,:,:);

% allmb_aroma_cor_srt = allmb_aroma_cor(indexpower,indexpower,:,:);
% allmb_aroma_sdsd_srt = allmb_aroma_sdsd(indexpower,indexpower,:,:);
% allmb_aroma_gg_srt = allmb_aroma_gg(indexpower,indexpower,:,:);
% allmb_aroma_sd_node_srt=allmb_aroma_sd_node(indexpower,:,:);


save([resultspath 'CorrSDpower264_9p.mat'],'allmb_cor_srt','allmb_sdsd_srt','allmb_gg_srt')
save([resultspath 'CorrSDpower264_9p_tf.mat'],'allmb_tf_cor_srt','allmb_tf_sdsd_srt','allmb_tf_gg_srt')
save([resultspath 'CorrSDpower264_9p_sm.mat'],'allmb_sm_cor_srt','allmb_sm_sdsd_srt','allmb_sm_gg_srt')
save([resultspath 'CorrSDpower264_9p_tf_sm.mat'],'allmb_tf_sm_cor_srt','allmb_tf_sm_sdsd_srt','allmb_tf_sm_gg_srt')
% save([resultspath 'CorrSDpower264_9p_aroma.mat'],'allmb_aroma_cor_srt','allmb_aroma_sdsd_srt','allmb_aroma_gg_srt')

%%% Create g-factor file for each subject:
allmb_gf_node = NaN(nNode,9,32);
allmb_tf_gf_node = NaN(nNode,9,32);
allmb_sm_gf_node = NaN(nNode,9,32);
allmb_tf_sm_gf_node = NaN(nNode,9,32);
% allmb_aroma_gf_node = NaN(nNode,9,32);

for j=1:9
    allmb_gf_node(:,j,:) =  allmb_sd_node(:,j,:)./allmb_sd_node(:,2,:);
    allmb_tf_gf_node(:,j,:) =  allmb_tf_sd_node(:,j,:)./allmb_tf_sd_node(:,2,:);
    allmb_sm_gf_node(:,j,:) =  allmb_sm_sd_node(:,j,:)./allmb_sm_sd_node(:,2,:);
    allmb_tf_sm_gf_node(:,j,:) =  allmb_tf_sm_sd_node(:,j,:)./allmb_tf_sm_sd_node(:,2,:);
%    allmb_aroma_gf_node(:,j,:) = allmb_aroma_sd_node(:,j,:)./allmb_aroma_sd_node(:,2,:);
end

gf_9p_25=quantile(allmb_gf_node,0.25,3);
gf_9p_50=quantile(allmb_gf_node,0.5,3);
gf_9p_75=quantile(allmb_gf_node,0.75,3);

writetable(array2table(gf_9p_50,'VariableNames',mblist),'~/risk_share/OptimalSMS_rsfMRI/Results/median_gfactor_264_9p.csv')

gf_9p = NaN(264,27);
gf_9p(:,1:3:27)=gf_9p_25;
gf_9p(:,2:3:27)=gf_9p_50;
gf_9p(:,3:3:27)=gf_9p_75;
writetable(array2table(gf_9p),'~/risk_share/OptimalSMS_rsfMRI/Results/allquartiles_gfactor_264_9p.csv')

temp1 = gf_9p_50 - gf_9p_25;
temp2 = gf_9p_75 - gf_9p_50;
[temp1(:,7),temp2(:,7)]
% some big asymmetries; use quartiles


% create files with subject level gfactors:
save([resultspath 'Gfactor_power264_9p__9p_tf.mat'],'allmb_gf_node','allmb_tf_gf_node','mblist');



% aside: create a nifti file with a separate volume for each power node:
power = niftiread('~/risk_share/OptimalSMS_rsfMRI/PowerAtlas/total_sphere.nii');  
powerarray = zeros(91*109*91,264);

for k=1:264
            n_ind=find(power==k);
            powerarray(find(power==k),k) = 1;
end
powerarray = reshape(powerarray,[91,109,91,264]);
tempinfo.ImageSize=[91 109 91 264];
tempinfo.Datatype='uint8';
tempinfo.BitsPerPixel=8;
niftiwrite(uint8(powerarray),'~/risk_share/OptimalSMS_rsfMRI/PowerAtlas/total_sphere_multivolume.nii',tempinfo);
