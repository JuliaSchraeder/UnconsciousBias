clear all                                                                   % clear workspace

addpath ('C:/Users/juhoffmann/Desktop/eeglab2022.1')                          % path for eeg lab
eeglab;  close;                                                            % open eeg lab and close GUI

% define datapath 
DataPath = ('C:/Users/juhoffmann/Desktop/EEG_BIDS/EEG_250Hz/Matlab');

% read in subjects 
subjects = dir(fullfile(DataPath, '*Remove_Bad_Intervals.mat')); 


idx = ismember({subjects.name}, {'sub-007_BackwardMask_Remove_Bad_Intervals.mat','sub-008_BackwardMask_Remove_Bad_Intervals.mat'...
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
,'sub-129_BackwardMask_Remove_Bad_Intervals.mat','sub-130_BackwardMask_Remove_Bad_Intervals.mat','sub-131_BackwardMask_Remove_Bad_Intervals.mat'...
,'sub-004_BackwardMask_Remove_Bad_Intervals.mat','sub-006_BackwardMask_Remove_Bad_Intervals.mat'...
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


MDD_subjects = subjects(idx_MDD);
MDD_list = {MDD_subjects.name}';

all_subjects = subjects(idx);                                          

%Nguyen, V. T., & Cunnington, R. (2014). The superior temporal sulcus and the N170 during face processing: single trial analysis of concurrent EEG–fMRI. NeuroImage, 86, 492-502.
channel1 = 15;                                                              % channel P7
channel2 = 16;                                                              % channel P8
channel3 = 59;                                                              % channel P07
channel4 = 60;                                                              % channel P08

% Range for N170 within ERP...
% Range von 150ms-200ms für N170 sind --> 0,25*ms=datapoints. 
range_min = 88;                                                             % 200ms Baseline+150ms = 350ms (87,5 datapoints)
range_max = 100;                                                            % 200ms Baseline+200ms = 400ms (100 datapoints)



%% Run the script for all patients 


N170 = {};

for i = 1:length(all_subjects)
    FileName = all_subjects(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath);                                            % read in EEG data
    subName = FileName(1:7);    % get SubName
    
    if ismember(FileName, MDD_list)
        Group = "MDD";
    else 
        Group = "HC";
    end 

    Erp = BackwardMask_getERP(EEG,channel1,channel2,channel3,channel4,range_min,range_max);

    % To plot ERPs
    erp_N170_neutral_strong(i,:) = Erp.n_strong;
    erp_N170_happy_strong(i,:) = Erp.h_strong;
    erp_N170_sad_strong(i,:) = Erp.s_strong;
    
    erp_N170_neutral_weak(i,:) = Erp.n_weak;
    erp_N170_happy_weak(i,:) = Erp.h_weak;
    erp_N170_sad_weak(i,:) = Erp.s_weak;
    
    erp_N170_happy(i,:) = Erp.happy;
    erp_N170_neutral(i,:) = Erp.neutral;
    erp_N170_sad(i,:) = Erp.sad;
    
    
    erp_N170_weak(i,:) = Erp.weak;
    erp_N170_strong(i,:) = Erp.strong;

    erp_times = Erp.times;


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

    N170{i,10} = Erp.happy_N170;
    N170{i,11} = Erp.sad_N170;
    N170{i,12} = Erp.neutral_N170;
    N170{i,13} = Group;

end       




%% P100

% Range von 80ms-130ms für P100 sind --> 0,25*ms=datapoints. 
range_min = 70;                                                             % 200ms Baseline+80ms=280ms(70 datapoints)
range_max = 83;                                                             % 200+130=330(82,5 datapoints)


%% Run the script for all patients 

P100 = {};


for i = 1:length(all_subjects)
    FileName = all_subjects(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath);                                            % read in EEG data
    subName = FileName(1:7);                                                % get SubName
    
    if ismember(FileName, MDD_list)
        Group = "MDD";
    else 
        Group = "HC";
    end 

    Erp_P100 = BackwardMask_getERP_P100(EEG,channel1,channel2,channel3,channel4,range_min,range_max);

    % To plot ERPs
    erp_P100_neutral_strong(i,:) = Erp_P100.n_strong;
    erp_P100_happy_strong(i,:) = Erp_P100.h_strong;
    erp_P100_sad_strong(i,:) = Erp_P100.s_strong;
    
    erp_P100_neutral_weak(i,:) = Erp_P100.n_weak;
    erp_P100_happy_weak(i,:) = Erp_P100.h_weak;
    erp_P100_sad_weak(i,:) = Erp_P100.s_weak;
    
    erp_P100_happy(i,:) = Erp_P100.happy;
    erp_P100_neutral(i,:) = Erp_P100.neutral;
    erp_P100_sad(i,:) = Erp_P100.sad;
    
    
    erp_P100_weak(i,:) = Erp_P100.weak;
    erp_P100_strong(i,:) = Erp_P100.strong;

    erp_P100_times = Erp_P100.times;


    % For EEG-fMRI Analysis
    P100{i,1} = subName;

    P100{i,2} = Erp_P100.h_strong_P100;
    P100{i,3} = Erp_P100.n_strong_P100;
    P100{i,4} = Erp_P100.s_strong_P100;

    P100{i,5} = Erp_P100.h_weak_P100;
    P100{i,6} = Erp_P100.n_weak_P100;
    P100{i,7} = Erp_P100.s_weak_P100;

    P100{i,8} = Erp_P100.weak_P100;
    P100{i,9} = Erp_P100.strong_P100;

    P100{i,10} = Erp_P100.happy_P100;
    P100{i,11} = Erp_P100.sad_P100;
    P100{i,12} = Erp_P100.neutral_P100;
    P100{i,13} = Group;
end        
 

%% save

header = {'subName',...
    'h_strong','n_strong','s_strong',...
    'h_weak','n_weak','s_weak',...
    'weak','strong','happy','sad','neutral','group'};
P100df = [header; P100]; 
N170df = [header; N170];
writecell(P100df,'C:/Users/juhoffmann/Desktop/Git/UnconsciousBias/data/erp_P100_All.xlsx');
writecell(N170df,'C:/Users/juhoffmann/Desktop/Git/UnconsciousBias/data/erp_N170_All.xlsx');




%%
figure,
% sd for strong
y = mean(erp_N170_strong,1);
x = erp_times;
std_dev = std(erp_N170_strong);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0 0.4470 0.7410],'LineStyle',"none",'FaceAlpha',.3); hold on;
plot(x, mean(erp_N170_strong,1),'Color',[0 0.4470 0.7410],'LineWidth',1.5,'LineStyle','--'); 

%sd for weak
y = mean(erp_N170_weak,1);
x = erp_times;
std_dev = std(erp_N170_weak);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'magenta','LineStyle',"none",'FaceAlpha',.3); 
plot(x,mean(erp_N170_weak,1),'Color','magenta','LineWidth',1.5); hold off;

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('unconscious sd','unconscious','conscious sd','conscious','Location','northeast');
legend('boxoff'); title(['unconscious vs. conscious trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
ylim([-4 15])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'unconscious vs conscious trials.tiff');








%%
figure,
% sd for happy
y = mean(erp_N170_happy,1);
x = erp_times;
std_dev = std(erp_N170_happy);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.4660 0.6740 0.1880],'LineStyle',"none",'FaceAlpha',.5); hold on;
plot(x, mean(erp_N170_happy,1),'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5); 

%sd for sad
y = mean(erp_N170_sad,1);
x = erp_times;
std_dev = std(erp_N170_sad);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.9290 0.6940 0.1250],'LineStyle',"none",'FaceAlpha',.3); 
plot(x,mean(erp_N170_sad,1),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.5); 


%sd for neutral
y = mean(erp_N170_neutral,1);
x = erp_times;
std_dev = std(erp_N170_neutral);
curve1 = y + std_dev;
curve2 = y - std_dev;
x2 = [x, fliplr(x)];
inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, [0.3010 0.7450 0.9330],'LineStyle',"none",'FaceAlpha',.3); 
plot(x,mean(erp_N170_neutral,1),'Color',[0.3010 0.7450 0.9330],'LineWidth',1.5); hold off;

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
legend('happy sd','happy','sad sd','sad','neutral sd','neutral','Location','northeast');
legend('boxoff'); title(['unconscious vs. conscious trials ', 'P7/P07+P8/P08']);  
xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
ylim([-4 15])
set(gcf,'Position',[15 15 750 600])
print(gcf, '-dtiff', 'emotions.tiff');

















%%




% figure,
% plot(erp_times,mean(erp_strong,1),'Color','#0072BD','LineWidth',1.5,'LineStyle','--'); hold on; 
% plot(erp_times,mean(erp_weak,1),'Color','#0072BD','LineWidth',1.5); hold off;
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('unconscious','conscious','Location','northeast');
% legend('boxoff'); title(['unconscious vs. conscious trials ', 'P7/P07+P8/P08']);  
% xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
% ylim([-3 10])
% set(gcf,'Position',[15 15 750 600])
% print(gcf, '-dtiff', 'unconscious vs conscious trials.tiff');
% 
% 
% 
% 
% figure,
% plot(erp_times,mean(erp_happy,1),'Color','#77AC30','LineWidth',1.5,'LineStyle',':'); hold on; 
% plot(erp_times,mean(erp_neutral,1),'Color','#4DBEEE','LineWidth',1.5,'LineStyle','-.'); 
% plot(erp_times,mean(erp_sad,1),'Color','#EDB120','LineWidth',1.5); hold off;
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% legend('happy','neutral','sad','Location','northeast');
% legend('boxoff'); title(['unconscious vs. conscious trials ', 'P7/P07+P8/P08']);  
% xlim([-250 820]); xlabel('Amplitude (uV)', 'Position', [0 9],'FontSize', 8); ylabel('Time(ms)', 'Position',[0 0],'FontSize', 8); grid off; box off
% ylim([-3 10])
% set(gcf,'Position',[15 15 750 600])
% print(gcf, '-dtiff', 'emotions.tiff');



