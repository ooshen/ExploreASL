% This is a script created by Tamas Josza, for extracting ExploreASL CBF
% stats for his virtual brain model
%
% Suggestions from Henk Mutsaerts for simplification
% To be checked by Jan Petr for correctness PVC


clear all;
close all;
clc;

%% define files, folders, and kernel size
% TODO: please specify folders and Kernel size
wdir = './TestDataSet/';

% list of folders with ASL data
patient_folders = {'Sub-001/'};
ASL_data_folder = 'ASL_1/';

outfile = 'patient_stats.csv';

myKernel = [3,3,3];

%% initialise xASL
xASL_path = '/home/atom/Downloads/ExploreASL/';
addpath(xASL_path);

%fullfile('Development','dicomtools'), ...
subfolders_to_add = {'Functions', 'mex', 'Modules', ...
                        fullfile('Modules', 'SubModule_Structural'), ...
                        fullfile('Modules', 'SubModule_ASL'), ...
                        fullfile('Modules', 'SubModule_Population'), ...
                        'Development', 'External', ...
                        fullfile('External','isnear'), ...
                        fullfile('External','DCMTK'), ...
                        fullfile('External','ExploreQC'), ...
                        fullfile('External','SPMmodified'), ...
                        fullfile('External','SPMmodified','xASL'),...                            
                        fullfile('External','SPMmodified','toolbox','cat12'), ...
                        fullfile('External','SPMmodified','toolbox','LST'), ...
                        fullfile('External','SPMmodified','toolbox','OldNorm')};

for ii=1:length(subfolders_to_add)
    addpath(fullfile(xASL_path,subfolders_to_add{ii}));
end

%% compute statistical data for patients

% ID, age, gender, Vol_GM, Vol_WM, mCBF_GM, minCBF_GM, maxCBF_GM, mCBF_WM, minCBF_WM, maxCBF_WM
patient_data = zeros(length(patient_folders),11);
for ii = 1:length(patient_folders)
    
    load([wdir, patient_folders{1}, 'x.mat'])
    
    % determine voxel volume [ml]
    % TODO: please check voxel volume computation
    % (I do not know the difference between optimFWHM_Res_mm and optimFWHM)
    V_voxel = prod(x.S.optimFWHM_Res_mm)/1000;
    
    % read age, gender info
    % TODO: please implement reading patients' age and gender
    age = 70;
    gender = 0; % 0-female, 1-male
    
    % read images
    PVgm = xASL_io_Nifti2Im([wdir, patient_folders{ii}, ASL_data_folder, 'PVgm.nii.gz']);
    PVwm = xASL_io_Nifti2Im([wdir, patient_folders{ii}, ASL_data_folder, 'PVwm.nii.gz']);

    CBF = xASL_io_Nifti2Im([wdir, patient_folders{ii}, ASL_data_folder, 'CBF.nii.gz']);
    MaskVascular = xASL_io_Nifti2Im([wdir, patient_folders{ii}, ASL_data_folder, 'MaskVascular.nii.gz']);

    % check readings
%     [nx, ny, nz] = size(MaskVascular);
%
%     figure(1)
%     imshow(rot90(squeeze(PVgm(:,round(ny/2),:)),1))

    % compute perfusion stats [ml/min/100g]
    [meanCBF_GM_corr,   meanCBF_WM_corr]   = xASL_stat_ComputeMean(CBF,MaskVascular,[],2,PVgm,PVwm);

    [CBF_stack_corr,residual] = xASL_im_PVCkernel(CBF, cat(4,PVgm,PVwm), myKernel, 'gauss');
    CBF_GM_corr = CBF_stack_corr(:,:,:,1);
    CBF_WM_corr = CBF_stack_corr(:,:,:,2);

    GM_mask = PVgm>0.5 & MaskVascular==1;
    WM_mask = PVwm>0.5 & MaskVascular==1;

    CBF_GM_array = CBF_GM_corr(GM_mask);
    CBF_WM_array = CBF_WM_corr(WM_mask);
    
    % compute GM and WM volumes [ml]
    VGM = sum(V_voxel*PVgm(PVgm>0));
    VWM = sum(V_voxel*PVwm(PVwm>0));
    
    % create database
    patient_data(ii,:) = [ii, age, gender, VGM, VWM, ...
                          meanCBF_GM_corr, quantile(CBF_GM_array,[0.05, 0.95]), ...
                          meanCBF_WM_corr, quantile(CBF_WM_array,[0.05, 0.95])];
    
end

%% save database
patient_data = array2table(patient_data,'VariableNames',...
               {'ID','age','gender','Vol_GM','Vol_WM','mCBF_GM','p5CBF_GM','p95CBF_GM','mCBF_WM','p5CBF_WM','p95CBF_WM'});
writetable(patient_data,outfile);