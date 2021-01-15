function [] = SeedMaps_create(preprocList, roiName)
% BRisk
% Create seed maps
%       1) average correlation
%       2) Cohen's d
%
%
seedRegion=load_nii(strcat('~/risk_jhu_share/OptimalSMS_rsfMRI/Seeds/', roiName, '.nii.gz'));
nROI=1;

mask = load_nii('~/risk_jhu_share/OptimalSMS_rsfMRI/tissuepriors/MNI152_T1_2mm_brain_mask_ero.nii.gz');
inBrain = find(mask.img);

% check the size of the seed region:
% tabulate(seedRegion.img(:))

mblist={'rsfMRI_3_3mm_iso_MB1_AP','rsfMRI_2mm_iso_MB1_AP','rsfMRI_2mm_iso_MB2_AP','rsfMRI_2mm_iso_MB3_AP','rsfMRI_2mm_iso_MB4_AP',...
    'rsfMRI_2mm_iso_MB6_AP','rsfMRI_2mm_iso_MB8_AP', 'rsfMRI_2mm_iso_MB9_AP', 'rsfMRI_2mm_iso_MB12_AP'};

excludeScans = readmatrix('~/risk_jhu_share/OptimalSMS_rsfMRI/Results/exclude_subjects.csv');
subjList = excludeScans(:, 1);
includeScans = ~(excludeScans(:, 2:end));

nSubject = length(subjList);
nx=91;
ny=109;
nz=91;

% corrmatAll = nan(nx,ny,nz,nROI,9,nSubject);
corrmatAll = nan(nx,ny,nz, size(includeScans, 2), nSubject);

% note: sdmatAll made for just one roi...
% sdmatAll = nan(nx,ny,nz,9,nSubject);

for ipp = 1:length(preprocList)
    preproc = preprocList{ipp};
    template=load_nii(strcat('~/risk_jhu_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/102954/nii/rsfMRI_2mm_iso_MB1_AP/rsfMRI_2mm_iso_MB1_AP', preproc, '.nii.gz'));
    template.hdr.dime.glmax=1;
    template.hdr.dime.glmin=-1;
    template.hdr.dime.dim(5) = length(mblist);
    for j=1:nSubject
        subj = num2str(subjList(j, :));
        parfor iscan = 1:length(mblist)
            %         skip bad scans
            if(includeScans(j, iscan))
                try
                    im_file=['~/risk_jhu_share/OptimalSMS_rsfMRI/rsfmri_multiband_nii_dicom/' subj '/nii/' mblist{iscan} '/' mblist{iscan} preproc '.nii.gz'];
                    s2=load_nii(im_file);
                    [nx,ny,nz,nt] = size(s2.img);
                    roi_ts = zeros(nt,nROI);
                    for n=1:nROI %loop through nodes
                        n_ind = find(seedRegion.img);
                        myimg = reshape(s2.img,[nx*ny*nz,nt]);
                        roi_ts(:,n)=mean(myimg(n_ind,:),1);
                    end
                    %csvwrite(['~/Dropbox/OptimalSMS_rsfMRI/Data/' subj '/' subj '_thalamusROI_ts_' mblist{i} '.csv'],roi_ts);
                    corrmat = corr(roi_ts,reshape(permute(s2.img,[4,1,2,3]),[nt,nx*ny*nz]));
                    temp = reshape(corrmat,[nROI,nx,ny,nz]);
                    temp = permute(temp,[2,3,4,1]);
                    corrmatAll(:, :, :,iscan,j) = temp;
                    
                    fprintf(['Finished ' mblist{iscan} '\n*********\n'])
                catch
                    warning(['No ' mblist{iscan} preproc ' for participant ' subj])
                end
            end
        end
        fprintf(['Finished ' subj '\n*********\n'])
    end
    %193728 and 195590 missing MB9
    corrmatAll(repmat(mask.img,[1,1,1,9,32])==0) = NaN;
    
    % Optional: save subject levels maps
    %save(fullfile('~/risk_jhu_share/OptimalSMS_rsfMRI/Seeds', strcat(roiName, preproc, '_corrmatAll_excludedScans_mask_ero.mat')),'corrmatAll')
    
    % for one roi:
    meancormat = squeeze(nanmean(corrmatAll(:,:,:,:,:),5));
    sum(isnan(meancormat(:)))
    
    meancormat(repmat(mask.img,[1,1,1,9])==0) = NaN;
    sum(isnan(meancormat(:)))
    
    template.img = meancormat;
    
    save_nii(template,fullfile('~/risk_jhu_share/OptimalSMS_rsfMRI/Seeds', ...
        strcat(roiName, preproc, '_excludedScans_mask_ero_fc.nii.gz')))
    
    % estimate t-statistics:
    zcorrmatAll = atanh(corrmatAll);
    
    %Tabulate the number of voxels that have data for each acquisition.
    %The rows correspond to the number of subjects that have data
    %The first row equals 0 corresponds to mask=0, so most voxels have no
    %data
    %The last row is the most important: This count should be
    %compared to 194369. 
    % The maximum value (e.g., 31) correspond to the number of subects
    % with usable data for that acquisition.
    % The 2 mm acqusitions should have similar
    % numbers, however, the SB 3.3 mm tends to be a little lower due to
    % interpolation and registration to MNI 2 mm
    for t=1:9
        fprintf(['Counts for ',mblist{t}])
        fprintf('\n')
        mbcounts = sum(~isnan(corrmatAll(:,:,:,t,:)),5);
        tabulate(mbcounts(:))
        fprintf('***********************')
        fprintf('\n')
    end
    
    sdmap = nanstd(zcorrmatAll, 0, 5);
    template.img = sdmap;
    save_nii(template, fullfile('~/risk_jhu_share/OptimalSMS_rsfMRI/Seeds', ...
        strcat(roiName, preproc, '_excludedScans_mask_ero_sd.nii.gz')))
    
    % already masked because inputs are masked:
    dmap = nanmean(zcorrmatAll, 5)./sdmap;
    % sum(~isnan(dmap(:)))
    
    % na voxels with few counts
    % Currently, we use the masked image without this additional control
    % counts = sum(~isnan(corrmatAll),5);
    % dmap_2 = dmap;
    % dmap_2(counts<2) = NaN;
    
    % (sum(~isnan(dmap(:))) - sum(~isnan(dmap_2(:))))/9
    % lose ~47626 with <10, 33831 with <5, ~17000 with <2
    
    template.img = dmap;
    save_nii(template, fullfile('~/risk_jhu_share/OptimalSMS_rsfMRI/Seeds', ...
        strcat(roiName, preproc, '_excludedScans_mask_ero_cohensD.nii.gz')))
    
    
    
end
