clear all                                                                   % clear workspace

addpath ('C:/Users/juhoffmann/Desktop/eeglab2022.1')                          % path for eeg lab
eeglab;                                                              % open eeg lab and close GUI

% define datapath 
DataPath = ('C:/Users/juhoffmann/Desktop/EEG_BIDS/EEG_250Hz/Matlab');

% read in subjects 
subjects = dir(fullfile(DataPath, '*Remove_Bad_Intervals.mat')); 


% checked 11.04.2023
% Take MDD participants (exclusion: sub-048, sub-060, sub-061, sub-067, sub-087, sub-092, sub-095, sub-121)
idx_MDD = ismember({subjects.name}, {'sub-007_BackwardMask_Remove_Bad_Intervals.mat','sub-008_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-009_BackwardMask_Remove_Bad_Intervals.mat','sub-012_BackwardMask_Remove_Bad_Intervals.mat','sub-020_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-026_BackwardMask_Remove_Bad_Intervals.mat','sub-035_BackwardMask_Remove_Bad_Intervals.mat','sub-036_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-037_BackwardMask_Remove_Bad_Intervals.mat','sub-038_BackwardMask_Remove_Bad_Intervals.mat','sub-039_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-040_BackwardMask_Remove_Bad_Intervals.mat','sub-042_BackwardMask_Remove_Bad_Intervals.mat','sub-044_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-049_BackwardMask_Remove_Bad_Intervals.mat','sub-064_BackwardMask_Remove_Bad_Intervals.mat','sub-065_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-066_BackwardMask_Remove_Bad_Intervals.mat','sub-072_BackwardMask_Remove_Bad_Intervals.mat','sub-075_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-076_BackwardMask_Remove_Bad_Intervals.mat','sub-077_BackwardMask_Remove_Bad_Intervals.mat','sub-081_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-082_BackwardMask_Remove_Bad_Intervals.mat','sub-084_BackwardMask_Remove_Bad_Intervals.mat','sub-094_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-097_BackwardMask_Remove_Bad_Intervals.mat','sub-098_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-099_BackwardMask_Remove_Bad_Intervals.mat','sub-100_BackwardMask_Remove_Bad_Intervals.mat','sub-102_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-104_BackwardMask_Remove_Bad_Intervals.mat','sub-106_BackwardMask_Remove_Bad_Intervals.mat','sub-107_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-108_BackwardMask_Remove_Bad_Intervals.mat','sub-109_BackwardMask_Remove_Bad_Intervals.mat','sub-110_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-111_BackwardMask_Remove_Bad_Intervals.mat','sub-112_BackwardMask_Remove_Bad_Intervals.mat','sub-113_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-114_BackwardMask_Remove_Bad_Intervals.mat','sub-115_BackwardMask_Remove_Bad_Intervals.mat','sub-116_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-117_BackwardMask_Remove_Bad_Intervals.mat','sub-118_BackwardMask_Remove_Bad_Intervals.mat','sub-120_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-122_BackwardMask_Remove_Bad_Intervals.mat','sub-124_BackwardMask_Remove_Bad_Intervals.mat','sub-128_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-129_BackwardMask_Remove_Bad_Intervals.mat','sub-130_BackwardMask_Remove_Bad_Intervals.mat','sub-131_BackwardMask_Remove_Bad_Intervals.mat'});

mdd_subjects = subjects(idx_MDD);                                           % take only MDD subjects


% checked 11.04.2023
% Take HC participants (exclusion: sub-013, sub-016, sub-023, sub-028, sub-029, sub-030, sub-050, sub-055, sub-058, sub-063, sub-080, sub-119)
idx_HC = ismember({subjects.name}, {'sub-004_BackwardMask_Remove_Bad_Intervals.mat','sub-006_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-010_BackwardMask_Remove_Bad_Intervals.mat','sub-011_BackwardMask_Remove_Bad_Intervals.mat','sub-014_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-015_BackwardMask_Remove_Bad_Intervals.mat','sub-017_BackwardMask_Remove_Bad_Intervals.mat','sub-018_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-019_BackwardMask_Remove_Bad_Intervals.mat','sub-021_BackwardMask_Remove_Bad_Intervals.mat','sub-022_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-024_BackwardMask_Remove_Bad_Intervals.mat','sub-025_BackwardMask_Remove_Bad_Intervals.mat','sub-027_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-031_BackwardMask_Remove_Bad_Intervals.mat','sub-032_BackwardMask_Remove_Bad_Intervals.mat','sub-033_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-034_BackwardMask_Remove_Bad_Intervals.mat','sub-041_BackwardMask_Remove_Bad_Intervals.mat','sub-043_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-045_BackwardMask_Remove_Bad_Intervals.mat','sub-046_BackwardMask_Remove_Bad_Intervals.mat','sub-047_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-051_BackwardMask_Remove_Bad_Intervals.mat','sub-052_BackwardMask_Remove_Bad_Intervals.mat','sub-053_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-054_BackwardMask_Remove_Bad_Intervals.mat','sub-056_BackwardMask_Remove_Bad_Intervals.mat','sub-057_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-059_BackwardMask_Remove_Bad_Intervals.mat','sub-062_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-068_BackwardMask_Remove_Bad_Intervals.mat','sub-069_BackwardMask_Remove_Bad_Intervals.mat','sub-070_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-071_BackwardMask_Remove_Bad_Intervals.mat','sub-073_BackwardMask_Remove_Bad_Intervals.mat','sub-074_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-078_BackwardMask_Remove_Bad_Intervals.mat','sub-079_BackwardMask_Remove_Bad_Intervals.mat','sub-083_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-085_BackwardMask_Remove_Bad_Intervals.mat','sub-086_BackwardMask_Remove_Bad_Intervals.mat','sub-088_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-089_BackwardMask_Remove_Bad_Intervals.mat','sub-090_BackwardMask_Remove_Bad_Intervals.mat','sub-091_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-093_BackwardMask_Remove_Bad_Intervals.mat','sub-096_BackwardMask_Remove_Bad_Intervals.mat','sub-101_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-103_BackwardMask_Remove_Bad_Intervals.mat','sub-105_BackwardMask_Remove_Bad_Intervals.mat','sub-123_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-125_BackwardMask_Remove_Bad_Intervals.mat','sub-126_BackwardMask_Remove_Bad_Intervals.mat','sub-127_BackwardMask_Remove_Bad_Intervals.mat'});

healthy_subjects = subjects(idx_HC);                                        % take only HC subjects


%Nguyen, V. T., & Cunnington, R. (2014). The superior temporal sulcus and the N170 during face processing: single trial analysis of concurrent EEG–fMRI. NeuroImage, 86, 492-502.
channel1 = 15;                                                              % channel P7
channel2 = 16;                                                              % channel P8
channel3 = 59;                                                              % channel P07
channel4 = 60;                                                              % channel P08

% Range for N170 within ERP...
% Range von 150ms-200ms für N170 sind --> 0,25*ms=datapoints. 
range_min = 88;                                                             % 200ms Baseline+150ms = 350ms (87,5 datapoints)
range_max = 100;                                                            % 200ms Baseline+200ms = 400ms (100 datapoints)



%% MDD: Run the script for MDD patients 

N170 = {};
dataframe = {};
dataframe_conscious_effect = {};

for i = 1:length(mdd_subjects)
    FileName = mdd_subjects(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath);                                            % read in EEG data
    subName = FileName(1:7);                                                % get SubName
    
    Erp = BackwardMask_getERP(EEG,channel1,channel2,channel3,channel4,range_min,range_max);

    % To plot ERPs
    erp_neutral_strong(i,:) = Erp.n_strong;
    erp_happy_strong(i,:) = Erp.h_strong;
    erp_sad_strong(i,:) = Erp.s_strong;
    
    erp_neutral_weak(i,:) = Erp.n_weak;
    erp_happy_weak(i,:) = Erp.h_weak;
    erp_sad_weak(i,:) = Erp.s_weak;

    erp_weak(i,:) = Erp.weak;
    erp_strong(i,:) = Erp.strong;

    erp_times = Erp.times;

        % Effects
        
    erp_cons_effect(i,:) = Erp.cons_effect;
    erp_cons_effect_happy(i,:) = Erp.cons_effect_happy;
    erp_cons_effect_neutral(i,:) = Erp.cons_effect_neutral;
    erp_cons_effect_sad(i,:) = Erp.cons_effect_sad;

    erp_emotion_effect_happy_conscious(i,:) = Erp.emotion_effect_happy_conscious;
    erp_emotion_effect_sad_conscious(i,:) = Erp.emotion_effect_sad_conscious;
    erp_emotion_effect_happy_unconscious(i,:) = Erp.emotion_effect_happy_unconscious;
    erp_emotion_effect_sad_unconscious(i,:) = Erp.emotion_effect_sad_unconscious;

    
    % For EEG-fMRI Analysis
    N170{i,1} = subName;

    N170{i,2} = Erp.h_strong_N170;
    N170{i,3} = Erp.n_strong_N170;
    N170{i,4} = Erp.s_strong_N170;

    N170{i,5} = Erp.h_weak_N170;
    N170{i,6} = Erp.n_weak_N170;
    N170{i,7} = Erp.s_weak_N170;

    N170{i,8} = Erp.weak_N170;
    N170{i,9} = Erp.strong_N170;

 
    N170{i,10} = Erp.h_h_weak_N170; 
    N170{i,11} = Erp.h_n_weak_N170; 
    N170{i,12} = Erp.h_s_weak_N170;
  
    N170{i,13} = Erp.s_h_weak_N170; 
    N170{i,14} = Erp.s_n_weak_N170; 
    N170{i,15} = Erp.s_s_weak_N170; 
  
    N170{i,16} = Erp.n_h_weak_N170; 
    N170{i,17} = Erp.n_n_weak_N170; 
    N170{i,18} = Erp.n_s_weak_N170; 
      
    N170{i,19} = Erp.h_h_strong_N170; 
    N170{i,20} = Erp.h_n_strong_N170; 
    N170{i,21} = Erp.h_s_strong_N170; 
      
    N170{i,22} = Erp.s_h_strong_N170; 
    N170{i,23} = Erp.s_n_strong_N170; 
    N170{i,24} = Erp.s_s_strong_N170; 
      
    N170{i,25} = Erp.n_h_strong_N170; 
    N170{i,26} = Erp.n_n_strong_N170; 
    N170{i,27} = Erp.n_s_strong_N170; 
    
        % Effects
        
    N170{i,28} = Erp.cons_effect_N170; 
    N170{i,29} = Erp.cons_effect_happy_N170; 
    N170{i,30} = Erp.cons_effect_neutral_N170; 
    N170{i,31} = Erp.cons_effect_sad_N170; 
    
    N170{i,32} = Erp.emotion_effect_happy_conscious_N170; 
    N170{i,33} = Erp.emotion_effect_sad_conscious_N170; 
    N170{i,34} = Erp.emotion_effect_happy_unconscious_N170; 
    N170{i,35} = Erp.emotion_effect_sad_unconscious_N170;    
    
    
    %% Längsdaten

    j=1+(6*(i-1));    
    
    k=j+1;
    l=k+1;
    m=l+1;
    n=m+1;
    o=n+1;
    
    dataframe{j,1}= subName;
    dataframe{k,1}= subName;
    dataframe{l,1}= subName;
    dataframe{m,1}= subName;
    dataframe{n,1}= subName;
    dataframe{o,1}= subName;
    
    dataframe{j,2}= 'bewusst';
    dataframe{k,2}= 'bewusst';
    dataframe{l,2}= 'bewusst';
    dataframe{m,2}= 'unbewusst';
    dataframe{n,2}= 'unbewusst';
    dataframe{o,2}= 'unbewusst';
    
    dataframe{j,3}= 'happy';
    dataframe{k,3}= 'sad';
    dataframe{l,3}= 'neutral';
    dataframe{m,3}= 'happy';
    dataframe{n,3}= 'sad';
    dataframe{o,3}= 'neutral';
    
    dataframe{j,5}= 'happy_bewusst';
    dataframe{j,4}= Erp.h_weak_N170;

    dataframe{k,5}= 'sad_bewusst';
    dataframe{k,4}= Erp.s_weak_N170;
    
    dataframe{l,5}= 'neutral_bewusst';
    dataframe{l,4}= Erp.n_weak_N170;
    
    dataframe{m,5}= 'happy_unbewusst';
    dataframe{m,4}= Erp.h_strong_N170;
    
    dataframe{n,5}='sad_unbewusst';
    dataframe{n,4}= Erp.s_strong_N170;
    
    dataframe{o,5}='neutral_unbewusst';
    dataframe{o,4}= Erp.n_strong_N170; 
    
    dataframe{j,6}= 'MDD';
    dataframe{k,6}= 'MDD';
    dataframe{l,6}= 'MDD';
    dataframe{m,6}= 'MDD';
    dataframe{n,6}= 'MDD';
    dataframe{o,6}= 'MDD';
    
    %% Effects

    j=1+(3*(i-1));    
    
    k=j+1;
    l=k+1;

    
    dataframe_conscious_effect{j,1}= subName;
    dataframe_conscious_effect{k,1}= subName;
    dataframe_conscious_effect{l,1}= subName;
    
    dataframe_conscious_effect{j,2}= 'happy';
    dataframe_conscious_effect{k,2}= 'sad';
    dataframe_conscious_effect{l,2}= 'neutral';

    dataframe_conscious_effect{j,3}= Erp.cons_effect_happy_N170;
    dataframe_conscious_effect{k,3}= Erp.cons_effect_sad_N170;
    dataframe_conscious_effect{l,3}= Erp.cons_effect_neutral_N170;
    
    dataframe_conscious_effect{j,4}= 'MDD';
    dataframe_conscious_effect{k,4}= 'MDD';
    dataframe_conscious_effect{l,4}= 'MDD';

    

end
% Längs

header = {'subName','emotion','mask','N170','condition','group'};
dataframe = [header; dataframe]; 
writecell(dataframe,'.../erp_N170_MDD_Längs.xlsx');

% Längs Conscious effect
header = {'subName','emotion','N170','group'};
dataframe_conscious_effect = [header; dataframe_conscious_effect]; 
writecell(dataframe_conscious_effect,'.../erp_N170_MDD_Längs_conscious_effect.xlsx');


% Quer
header = {'subName',...
    'h_strong_N170','n_strong_N170','s_strong_N170',...
    'h_weak_N170','n_weak_N170','s_weak_N170',...
    'weak_N170','strong_N170',...
    'h_h_weak_N170','h_n_weak_N170','h_s_weak_N170',...
    's_h_weak_N170','s_n_weak_N170','s_s_weak_N170',...
    'n_h_weak_N170','n_n_weak_N170','n_s_weak_N170',...
    'h_h_strong_N170','h_n_strong_N170','h_s_strong_N170',...
    's_h_strong_N170','s_n_strong_N170','s_s_strong_N170',...
    'n_h_strong_N170','n_n_strong_N170','n_s_strong_N170',...
    'cons_effect','cons_effect_happy','cons_effect_sad','cons_effect_neutral',...
    'emotion_effect_happy_conscious','emotion_effect_sad_conscious',...
    'emotion_effect_happy_unconscious','emotion_effect_sad_unconscious'};

N170 = [header; N170];
writecell(N170, '.../erp_N170_MDD_Quer.xlsx', 'Sheet', 1, 'Range', 'D1')


%% HC: Run the script for HC participants

N170 = {};
dataframe = {};
dataframe_conscious_effect = {};

for i = 1:length(healthy_subjects)
    FileName = healthy_subjects(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath);                                            % read in EEG data
    subName = FileName(1:7);                                                % get SubName
    
    Erp = BackwardMask_getERP(EEG,channel1,channel2,channel3,channel4,range_min,range_max);

    % To plot ERPs
    HC_erp_neutral_strong(i,:) = Erp.n_strong;
    HC_erp_happy_strong(i,:) = Erp.h_strong;
    HC_erp_sad_strong(i,:) = Erp.s_strong;
    
    HC_erp_neutral_weak(i,:) = Erp.n_weak;
    HC_erp_happy_weak(i,:) = Erp.h_weak;
    HC_erp_sad_weak(i,:) = Erp.s_weak;

    HC_erp_weak(i,:) = Erp.weak;
    HC_erp_strong(i,:) = Erp.strong;

    HC_erp_times = Erp.times;

            % Effects
        
    HC_erp_cons_effect(i,:) = Erp.cons_effect;
    HC_erp_cons_effect_happy(i,:) = Erp.cons_effect_happy;
    HC_erp_cons_effect_neutral(i,:) = Erp.cons_effect_neutral;
    HC_erp_cons_effect_sad(i,:) = Erp.cons_effect_sad;

    HC_erp_emotion_effect_happy_conscious(i,:) = Erp.emotion_effect_happy_conscious;
    HC_erp_emotion_effect_sad_conscious(i,:) = Erp.emotion_effect_sad_conscious;
    HC_erp_emotion_effect_happy_unconscious(i,:) = Erp.emotion_effect_happy_unconscious;
    HC_erp_emotion_effect_sad_unconscious(i,:) = Erp.emotion_effect_sad_unconscious;
    

    % For EEG-fMRI Analysis
    N170{i,1} = subName;

    N170{i,2} = Erp.h_strong_N170;
    N170{i,3} = Erp.n_strong_N170;
    N170{i,4} = Erp.s_strong_N170;

    N170{i,5} = Erp.h_weak_N170;
    N170{i,6} = Erp.n_weak_N170;
    N170{i,7} = Erp.s_weak_N170;

    N170{i,8} = Erp.weak_N170;
    N170{i,9} = Erp.strong_N170;

 
    N170{i,10} = Erp.h_h_weak_N170; 
    N170{i,11} = Erp.h_n_weak_N170; 
    N170{i,12} = Erp.h_s_weak_N170;
  
    N170{i,13} = Erp.s_h_weak_N170; 
    N170{i,14} = Erp.s_n_weak_N170; 
    N170{i,15} = Erp.s_s_weak_N170; 
  
    N170{i,16} = Erp.n_h_weak_N170; 
    N170{i,17} = Erp.n_n_weak_N170; 
    N170{i,18} = Erp.n_s_weak_N170; 
      
    N170{i,19} = Erp.h_h_strong_N170; 
    N170{i,20} = Erp.h_n_strong_N170; 
    N170{i,21} = Erp.h_s_strong_N170; 
      
    N170{i,22} = Erp.s_h_strong_N170; 
    N170{i,23} = Erp.s_n_strong_N170; 
    N170{i,24} = Erp.s_s_strong_N170; 
      
    N170{i,25} = Erp.n_h_strong_N170; 
    N170{i,26} = Erp.n_n_strong_N170; 
    N170{i,27} = Erp.n_s_strong_N170; 
        
        % Effects
        
    N170{i,28} = Erp.cons_effect_N170; 
    N170{i,29} = Erp.cons_effect_happy_N170; 
    N170{i,30} = Erp.cons_effect_neutral_N170; 
    N170{i,31} = Erp.cons_effect_sad_N170; 
    
    N170{i,32} = Erp.emotion_effect_happy_conscious_N170; 
    N170{i,33} = Erp.emotion_effect_sad_conscious_N170; 
    N170{i,34} = Erp.emotion_effect_happy_unconscious_N170; 
    N170{i,35} = Erp.emotion_effect_sad_unconscious_N170;    
    
    
    %% Längsdaten

    j=1+(6*(i-1));    
    
    k=j+1;
    l=k+1;
    m=l+1;
    n=m+1;
    o=n+1;
    
    dataframe{j,1}= subName;
    dataframe{k,1}= subName;
    dataframe{l,1}= subName;
    dataframe{m,1}= subName;
    dataframe{n,1}= subName;
    dataframe{o,1}= subName;
    
    dataframe{j,2}= 'bewusst';
    dataframe{k,2}= 'bewusst';
    dataframe{l,2}= 'bewusst';
    dataframe{m,2}= 'unbewusst';
    dataframe{n,2}= 'unbewusst';
    dataframe{o,2}= 'unbewusst';
    
    dataframe{j,3}= 'happy';
    dataframe{k,3}= 'sad';
    dataframe{l,3}= 'neutral';
    dataframe{m,3}= 'happy';
    dataframe{n,3}= 'sad';
    dataframe{o,3}= 'neutral';
    
    dataframe{j,5}= 'happy_bewusst';
    dataframe{j,4}= Erp.h_weak_N170;

    dataframe{k,5}= 'sad_bewusst';
    dataframe{k,4}= Erp.s_weak_N170;
    
    dataframe{l,5}= 'neutral_bewusst';
    dataframe{l,4}= Erp.n_weak_N170;
    
    dataframe{m,5}= 'happy_unbewusst';
    dataframe{m,4}= Erp.h_strong_N170;
    
    dataframe{n,5}='sad_unbewusst';
    dataframe{n,4}= Erp.s_strong_N170;
    
    dataframe{o,5}='neutral_unbewusst';
    dataframe{o,4}= Erp.n_strong_N170; 
    
    dataframe{j,6}= 'HC';
    dataframe{k,6}= 'HC';
    dataframe{l,6}= 'HC';
    dataframe{m,6}= 'HC';
    dataframe{n,6}= 'HC';
    dataframe{o,6}= 'HC';
    
    %% Effects

    j=1+(3*(i-1));    
    
    k=j+1;
    l=k+1;

    
    dataframe_conscious_effect{j,1}= subName;
    dataframe_conscious_effect{k,1}= subName;
    dataframe_conscious_effect{l,1}= subName;
    
    dataframe_conscious_effect{j,2}= 'happy';
    dataframe_conscious_effect{k,2}= 'sad';
    dataframe_conscious_effect{l,2}= 'neutral';

    dataframe_conscious_effect{j,3}= Erp.cons_effect_happy_N170;
    dataframe_conscious_effect{k,3}= Erp.cons_effect_sad_N170;
    dataframe_conscious_effect{l,3}= Erp.cons_effect_neutral_N170;
    
    dataframe_conscious_effect{j,4}= 'HC';
    dataframe_conscious_effect{k,4}= 'HC';
    dataframe_conscious_effect{l,4}= 'HC';
    
    
end
% Längs

header = {'subName','emotion','mask','N170','condition','group'};
dataframe = [header; dataframe]; 
writecell(dataframe,'.../erp_N170_HC_Längs.xlsx');

% Längs Conscious effect
header = {'subName','emotion','N170','group'};
dataframe_conscious_effect = [header; dataframe_conscious_effect]; 
writecell(dataframe_conscious_effect,'.../erp_N170_HC_Längs_conscious_effect.xlsx');


% Quer
header = {'subName',...
    'h_strong_N170','n_strong_N170','s_strong_N170',...
    'h_weak_N170','n_weak_N170','s_weak_N170',...
    'weak_N170','strong_N170',...
    'h_h_weak_N170','h_n_weak_N170','h_s_weak_N170',...
    's_h_weak_N170','s_n_weak_N170','s_s_weak_N170',...
    'n_h_weak_N170','n_n_weak_N170','n_s_weak_N170',...
    'h_h_strong_N170','h_n_strong_N170','h_s_strong_N170',...
    's_h_strong_N170','s_n_strong_N170','s_s_strong_N170',...
    'n_h_strong_N170','n_n_strong_N170','n_s_strong_N170',...
    'cons_effect','cons_effect_happy','cons_effect_sad','cons_effect_neutral',...
    'emotion_effect_happy_conscious','emotion_effect_sad_conscious',...
    'emotion_effect_happy_unconscious','emotion_effect_sad_unconscious'};

N170 = [header; N170];
writecell(N170, '.../erp_N170_HC_Quer.xlsx', 'Sheet', 1, 'Range', 'D1')


%%

%% Plot %%

%% HC vs MDD

figure,
plot(erp_times,mean(erp_strong,1),'Color','#FF0000','LineWidth',1.5,'LineStyle','--'); hold on;
plot(erp_times,mean(HC_erp_strong,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); 
plot(erp_times,mean(erp_weak,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(HC_erp_weak,1),'Color','#0072BD','LineWidth',1.5); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('MDD unconscious','HC unconscious','MDD conscious','HC conscious','Location','northeast');
legend('boxoff'); title(['HC vs MDD: weakly vs strongly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
ylim([-3 10])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'HC vs MDD weakly vs strongly masked trials.tiff');



% only HC

figure,

plot(erp_times,mean(HC_erp_strong,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold on;
plot(erp_times,mean(HC_erp_weak,1),'Color','#0072BD','LineWidth',1.5); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('unconscious','conscious','Location','northeast');
legend('boxoff'); title(['HC: weakly vs strongly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
ylim([-3 10])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'HC weakly vs strongly masked trials.tiff');

% only MDD

figure,

plot(erp_times,mean(erp_strong,1),'Color','#FF0000','LineWidth',1.5,'LineStyle','--'); hold on;
plot(erp_times,mean(erp_weak,1),'Color','#FF0000','LineWidth',1.5); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('unconscious','conscious','Location','northeast');
legend('boxoff'); title(['MDD: weakly vs strongly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
ylim([-3 10])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'MDD weakly vs strongly masked trials.tiff');


%% Compare unconscious emotions

%MDD
figure,
plot(erp_times,mean(erp_happy_strong,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(erp_sad_strong,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(erp_neutral_strong,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northeast');
legend('boxoff'); title(['MDD: strongly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off 
ylim([-3 10])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'MDD strongly masked trials.tiff');




% HC
figure,
plot(erp_times,mean(HC_erp_happy_strong,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(HC_erp_sad_strong,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(HC_erp_neutral_strong,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northeast');
legend('boxoff'); title(['HC: strongly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off 
ylim([-3 10])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'HC strongly masked trials.tiff');


%% Compare concious emotions

%MDD
figure,
plot(erp_times,mean(erp_happy_weak,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(erp_sad_weak,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(erp_neutral_weak,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--');  hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northeast');
legend('boxoff'); title(['MDD: weakly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 6.5],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off 
ylim([-2 7])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'MDD weakly masked trials.tiff');


% HC
figure,
plot(erp_times,mean(HC_erp_happy_weak,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(HC_erp_sad_weak,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(HC_erp_neutral_weak,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northeast');
legend('boxoff'); title(['HC: weakly masked trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 6.5],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off  
set(gcf,'Position',[15 15 750 600])
ylim([-2 7])
print(gcf, '-dtiff', 'HC weakly masked trials.tiff');







%% Compare effects
% Consciousness

%MDD
figure,
plot(erp_times,mean(erp_cons_effect_happy,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(erp_cons_effect_sad,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(erp_cons_effect_neutral,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--');  hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northwest');
legend('boxoff'); title(['MDD: consciousness effect (conscious-unconscious) ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 2.5],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 -0.4],'FontSize', 8); grid off; box off 
ylim([-6 4])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'MDD consciousness effect (conscious-unconscious).tiff');


% HC
figure,
plot(erp_times,mean(HC_erp_cons_effect_happy,1),'Color','#77AC30','LineWidth',1.5); hold on;
plot(erp_times,mean(HC_erp_cons_effect_sad,1),'Color','#FF0000','LineWidth',1.5);
plot(erp_times,mean(HC_erp_cons_effect_neutral,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold off;
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy','sad','neutral','Location','northwest');
legend('boxoff'); title(['HC: consciousness effect (conscious-unconscious) ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 2.5],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 -0.4],'FontSize', 8); grid off; box off  
ylim([-6 4])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'HC consciousness effect (conscious-unconscious).tiff');
