function preproc_BackwardMask_Julia(subName)
addpath '/bif/storage/storage1/projects/emocon/Scripts'

%define Paths
processing_folder    = '/bif/storage/storage1/projects/emocon/';
study_folder   = '/bif/storage/storage1/projects/emocon/Data/BIDS';
script_folder  = '/bif/storage/storage1/projects/emocon/Scripts';


%Go to project folder
cd (processing_folder)
CWD            = pwd;

%create Outputpath for every subject 
if ~exist(fullfile(processing_folder, 'Preproc', 'BackwardMask', subName), 'dir')
mkdir(fullfile(processing_folder, 'Preproc', 'BackwardMask', subName));
end
preprocDir = fullfile(processing_folder, 'Preproc', 'BackwardMask', subName);


spm('Defaults','fMRI');
spm_jobman('initcfg');
clear matlabbatch

%select batch with fieldmapping, realignment, coregistration, normalisation
%and smooting
preprocfile    = fullfile(script_folder, 'Preprocessing.mat');
load(preprocfile);

%% Select epi files
epiPath = fullfile(study_folder,subName,'ses-001', 'func');
epiName = strcat(subName,'_ses-001_task-BackwardMask_run-001_bold.nii');
epiZipName = strcat(subName,'_ses-001_task-BackwardMask_run-001_bold.nii.gz');
gunzip(fullfile(epiPath, epiZipName));                                      %unzip Epi images

epiFileArrayAll = spm_select('expand', fullfile(epiPath,epiName));          %selects all EPI files
number_EpiFiles = length(epiFileArrayAll);                                  %get number of EPI files 
nessesary_epiFiles = cellstr(epiFileArrayAll(6:number_EpiFiles,:));         
first_epiFile = cellstr(epiFileArrayAll(6,:));                              %take only the very first EPI


%% Select anatomy file
anatomyPath = fullfile(study_folder,subName,'ses-001', 'anat');             
anatName = strcat(subName,'_ses-001_run-001_T1w.nii');
anatZipName = strcat(subName,'_ses-001_run-001_T1w.nii.gz');

gunzip(fullfile(anatomyPath,anatZipName));                                  %unzip the compressed nifti file
anatomyFile= spm_select('ExtFPList', anatomyPath,anatName);                 %select anatomy file

%% Select fieldmap files
greyfieldPath = fullfile(study_folder,subName,'ses-001', 'fmap');
magnName = strcat(subName,'_ses-001_magnitude1.nii');
magnitudeFile = spm_select('ExtFPList',greyfieldPath,magnName);
phaseName = strcat(subName,'_ses-001_phasediff.nii');
phaseFile = spm_select('ExtFPList',greyfieldPath,phaseName);


%% Define Batch Inputs
%Fieldmap
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase     = {phaseFile};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = {magnitudeFile};
matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi                   = first_epiFile; 
%Realignment
matlabbatch{2}.spm.spatial.realignunwarp.data.scans                               = cellstr(nessesary_epiFiles);
%Coregistration
matlabbatch{3}.spm.spatial.coreg.estwrite.source                                  = cellstr(anatomyFile);

%% Move files to Preprocessing Folder %%
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(1) = cfg_dep('Calculate VDM: Voxel displacement map (Subj 1, Session 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','vdmfile', '{}',{1}));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(2) = cfg_dep('Realign & Unwarp: Realignment Param File (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','rpfile'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(3) = cfg_dep('Realign & Unwarp: Unwarp Params File (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','dsfile'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(4) = cfg_dep('Realign & Unwarp: Unwarped Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','uwrfiles'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(5) = cfg_dep('Realign & Unwarp: Unwarped Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','meanuwr'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(6) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(7) = cfg_dep('Coregister: Estimate & Reslice: Resliced Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rfiles'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(8) = cfg_dep('Normalise: Estimate & Write: Deformation (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','def'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(9) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.files(10) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{5}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{6}.cfg_basicio.file_dir.file_ops.file_move.action.moveto = {preprocDir};


%save Batch in preproc Dir
cd(preprocDir)
preproc_filename = sprintf('preproc_%s.mat',subName);
save(preproc_filename, 'matlabbatch')

%run created Batch
spm_jobman('run', matlabbatch);

%move back to starting directory
cd(CWD);

