clear; close all;

%% Set up paths (USER DEPENDENT!)
procRoot='/ISIS/proc5';
addpath('/usr/local/spm/spm12'); %add SPM to path
addpath([procRoot '/software/GitHub/T2Star']); %add T2* functions to path
addpath([procRoot '/software/GitHub/utilities']); %add utility functions to path

opts.dicomDir='./dicom_eg/2_3D_ME_SPGR'; %dicom folder for multi-echo acquisition
opts.niftiDir='./nifti_R2s'; %folder to create nifti files
opts.mapDir='./maps_R2s'; %folder to create parameter maps
opts.threshold1=40; %exclude voxels where first echo is less intense than this (avoid fitting voxels in air)
opts.threshold2=40; %exclude echoes from fitting if less intense than this (avoid fitting noisy data points)

pipeline_R2s_convert(opts);
pipeline_R2s_create_map(opts);
