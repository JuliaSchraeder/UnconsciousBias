clear all;

DataPath = fullfile('/bif/storage/storage1/projects/emocon/FirstLevel_Primer/BackwardMask');
Subjects = dir(fullfile(DataPath, 's*')); % get subjects

spm_jobman('initcfg');
spm('defaults', 'FMRI');
global defaults;



 for i = 69:numel(Subjects)
     PatPath = fullfile(DataPath, Subjects(i).name);
     
      %Define contrasts
      matlabbatch{1}.spm.stats.con.spmmat = {fullfile(PatPath, 'SPM.mat')};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'happy_uncon'; 
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'sad_uncon'; 
matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 1];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'neutral_uncon';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [0 0 1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'happy_con';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.convec = [0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'sad_con';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.convec = [0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'neutral_con';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.convec = [0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'intro';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [0 0 0 0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';


matlabbatch{1}.spm.stats.con.delete = 0; % delete old contrasts yes = 1, delete no = 0

      spm_jobman('run', matlabbatch);
      clear matlabbatch;

 end