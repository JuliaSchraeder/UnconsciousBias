

%% Full Factorial Design HC/MDD x Primer Emotion x Conscious


%% Find FirstLevel Contrasts
clear                                                                                               % clear Workspace
Data_Path = fullfile('/bif/storage/storage1/projects/emocon/FirstLevel_Primer/BackwardMask/');          % set DataPath


% HC excluded participants: sub-005,sub-013, sub-016, sub-058, sub-063, sub-101, sub-105 (movement), sub-078 (checked 25.09.2023)
HC_Subjects = {'sub-004','sub-006','sub-010','sub-011','sub-014','sub-015','sub-017','sub-018','sub-019','sub-021','sub-022','sub-024','sub-025'...
,'sub-027','sub-028','sub-029','sub-030','sub-031','sub-032','sub-033','sub-034','sub-041','sub-043','sub-045','sub-046','sub-047','sub-050','sub-051'...
,'sub-052','sub-053','sub-054','sub-055','sub-056','sub-057','sub-059','sub-062','sub-068','sub-069','sub-070','sub-071','sub-073','sub-074','sub-079'...
,'sub-080','sub-083','sub-085','sub-086','sub-088','sub-089','sub-090','sub-091','sub-093','sub-096','sub-103','sub-119','sub-123'...
,'sub-125','sub-126','sub-127'};

% MDD excluded participants: sub-020, sub-042, sub-094 (movement), sub-060, sub-061(checked 25.09.2023)
MDD_Subjects = {'sub-007','sub-008','sub-009','sub-012','sub-026','sub-035','sub-036','sub-037','sub-038','sub-039','sub-040','sub-044','sub-048'...
,'sub-049','sub-064','sub-065','sub-066','sub-067','sub-072','sub-075','sub-076','sub-077','sub-081','sub-082','sub-084','sub-087','sub-092'...
,'sub-095','sub-097','sub-098','sub-099','sub-100','sub-102','sub-104','sub-106','sub-107','sub-108','sub-109','sub-110','sub-111','sub-112','sub-113','sub-114'...
,'sub-115','sub-116','sub-117','sub-118','sub-120','sub-121','sub-122','sub-124','sub-128','sub-129','sub-130','sub-131'};

% Get HC Contrasts
for i = 1:numel(HC_Subjects)                                                                           
    subName = HC_Subjects{i};  
    HC_happy_unconscious{i,:} = fullfile(Data_Path, subName,'con_0001.nii,1');
    HC_sad_unconscious{i,:} = fullfile(Data_Path, subName,'con_0002.nii,1');
    HC_neutral_unconscious{i,:} = fullfile(Data_Path, subName,'con_0003.nii,1');
    HC_happy_conscious{i,:} = fullfile(Data_Path, subName,'con_0004.nii,1');
    HC_sad_conscious{i,:} = fullfile(Data_Path, subName,'con_0005.nii,1');
    HC_neutral_conscious{i,:} = fullfile(Data_Path, subName,'con_0006.nii,1');
end


% Get Patient Contrasts
for i = 1:numel(MDD_Subjects)
    subName = MDD_Subjects{i};
    MDD_happy_unconscious{i,:} = fullfile(Data_Path, subName,'con_0001.nii,1');
    MDD_sad_unconscious{i,:} = fullfile(Data_Path, subName,'con_0002.nii,1');
    MDD_neutral_unconscious{i,:} = fullfile(Data_Path, subName,'con_0003.nii,1');
    MDD_happy_conscious{i,:} = fullfile(Data_Path, subName,'con_0004.nii,1');
    MDD_sad_conscious{i,:} = fullfile(Data_Path, subName,'con_0005.nii,1');
    MDD_neutral_conscious{i,:} = fullfile(Data_Path, subName,'con_0006.nii,1');
    
end



%% Define Design
matlabbatch{1}.spm.stats.factorial_design.dir = {'/bif/storage/storage1/projects/emocon/SecondLevel/FullFact/Primer_effect'};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'emotion';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 3;        % happy = 1, sad = 2, neutral = 3
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1;                                          
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;      % varianzengleichheit ja = 1 nein = 0
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'consciousness';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;        % conscious = 1, unconscious = 2
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;

matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'group';
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 2;        % HC = 1, Patient = 2  
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).ancova = 0;


%% HC %%

% Unconscious
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    2
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = HC_happy_unconscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    2
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = HC_sad_unconscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    2
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = HC_neutral_unconscious;

% Conscious
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [1
                                                                    1
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = HC_happy_conscious;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    1
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = HC_sad_conscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    1
                                                                    1];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = HC_neutral_conscious;



%% MDD %%

% Unconscious
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1
                                                                    2
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = MDD_happy_unconscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    2
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = MDD_sad_unconscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    2
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = MDD_neutral_unconscious;

% Conscious
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [1
                                                                    1
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans = MDD_happy_conscious;

matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [2
                                                                    1
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = MDD_sad_conscious;


matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [3
                                                                    1
                                                                    2];
matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = MDD_neutral_conscious;


%%
matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%Estimate created Batch
matlabbatch{2}.spm.stats.fmri_est.spmmat = {'/bif/storage/storage1/projects/emocon/SecondLevel/FullFact/Primer_effect/SPM.mat'};
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


% Run created Batch
spm_jobman('initcfg');
spm_jobman('run', matlabbatch);
