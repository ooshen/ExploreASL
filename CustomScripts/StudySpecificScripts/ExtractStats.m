% This is a script created by Tamas Josza, for extracting ExploreASL CBF
% stats for his virtual brain model
%
% Suggestions from Henk Mutsaerts for simplification
% To be checked by Jan Petr for correctness PVC


clear all;
close all;
clc;

%% 1. define files, folders, and kernel size
% TODO: please specify folders and Kernel size
wdir = './TestDataSet/'; % note that its better to copy this outside of ExploreASL, to separate code and data

% list of folders with ASL data
patient_folders = {'Sub-001/'};
ASL_data_folder = 'ASL_1/';

outfile = 'patient_stats.csv';

myKernel = [3,3,3];

%% 2. initialise xASL
xASL_path = '/home/atom/Downloads/ExploreASL/';
CurrentPath = pwd;
cd(xASL_path);
ExploreASL_Master('',0); % initialize ExploreASL, without processing data
cd(CurrentPath); % back to the initial path

%% 3. compute statistical data for patients

% ID, age, sex, Vol_GM, Vol_WM, mCBF_GM, minCBF_GM, maxCBF_GM, mCBF_WM, minCBF_WM, maxCBF_WM
patient_data = zeros(length(patient_folders),11);
for ii = 1:length(patient_folders)
    
    x = load(fullfile(wdir, patient_folders{1}, 'x.mat'),'-mat');
    x = x.x; % looks superfluous but this seems best practice in Matlab, as you specify where the loaded file info goes to
    
    % determine voxel volume [mL] % nagging I know but mL I think is SI
    % TODO: please check voxel volume computation
    % (I do not know the difference between optimFWHM_Res_mm and optimFWHM)
    Volume_voxel = prod(x.S.optimFWHM_Res_mm)/1000; % I like to spell out so there can be no confusion of what "V" means ;)
    % Jan can answer your question here, note that some is in mm, other in
    % voxels, we use this FWHM as an absolute effective spatial resolution,
    % but also sometimes as smoothing that needs to be applied to the
    % structural/T1w segmentations to get to the ASL resolution (Jan this
    % should be clear in our code)
    % Instead of this parameter, you can take the parameter straight from
    % the NIfTI header:
    PathCBF = fullfile(wdir, patient_folders{ii}, ASL_data_folder, 'CBF.nii');
    NIfTI = xASL_io_ReadNifti(PathCBF);
    Volume_voxel = prod(NIfTI.hdr.pixdim(2:4));
    
    % read age, sex info
    % TODO: please implement reading patients' age and sex -> H: Will do
    age = 70;
    sex = 2; % 1-male, 2-female % usually sex is referred to, gender means something different. 
    % Trying to avoid the number 0, counting from 1 (0 can have a specific
    % meaning, and isn't nice as covariate).
    
    % read images
    PVgm = xASL_io_Nifti2Im(fullfile(wdir, patient_folders{ii}, ASL_data_folder, 'PVgm.nii')); % use fullfile (OS independent)
    % note that the xASL_ prefixed scripts all deal transparently with % .nii(.gz)
    PVwm = xASL_io_Nifti2Im(fullfile(wdir, patient_folders{ii}, ASL_data_folder, 'PVwm.nii'));

    CBF = xASL_io_Nifti2Im(PathCBF);
    MaskVascular = xASL_io_Nifti2Im(fullfile(wdir, patient_folders{ii}, ASL_data_folder, 'MaskVascular.nii'));

    % check readings
%     [nx, ny, nz] = size(MaskVascular);
%
%     figure(1); imshow(rot90(squeeze(PVgm(:,round(ny/2),:)),1)) % I can
%     recommend the DIP (Delft Image Library), with e.g. dip_image()
%     instead of imshow(), very nice to scroll through slices/volumes etc,
%     and joking
%     we do have our own replacements of this image processing toolbox,
%     e.g. xASL_im_rotate()

    % compute perfusion stats [mL/min/100g]
    [meanCBF_GM_corr, meanCBF_WM_corr] = xASL_stat_ComputeMean(CBF, MaskVascular, [], 2, PVgm, PVwm);

    [CBF_stack_corr,residual] = xASL_im_PVCkernel(CBF, cat(4,PVgm,PVwm), myKernel, 'gauss');
    CBF_GM_corr = CBF_stack_corr(:,:,:,1);
    CBF_WM_corr = CBF_stack_corr(:,:,:,2);

    GM_mask = PVgm>0.5 & MaskVascular==1;
    WM_mask = PVwm>0.5 & MaskVascular==1;

    CBF_GM_array = CBF_GM_corr(GM_mask);
    CBF_WM_array = CBF_WM_corr(WM_mask);
    
    % compute GM and WM volumes [mL]
    VGM = sum(Volume_voxel*PVgm(PVgm>0)); % if you mask the PVC CBF above with e.g. pGM>0.5, then you want to do the same here for getting the volume right?
    VWM = sum(Volume_voxel*PVwm(PVwm>0));
    
    % create database
    patient_data(ii,:) = [ii, age, sex, VGM, VWM, ...
                          meanCBF_GM_corr, quantile(CBF_GM_array,[0.05, 0.95]), ...
                          meanCBF_WM_corr, quantile(CBF_WM_array,[0.05, 0.95])];
    
end

%% save database
patient_data = array2table(patient_data,'VariableNames',...
               {'ID','age','gender','Vol_GM','Vol_WM','mCBF_GM','p5CBF_GM','p95CBF_GM','mCBF_WM','p5CBF_WM','p95CBF_WM'});
writetable(patient_data,outfile);