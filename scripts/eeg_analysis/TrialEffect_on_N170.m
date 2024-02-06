% Add EEGLAB to MATLAB path (update this path to your EEGLAB installation)
addpath('C:/Users/juhoffmann/Desktop/eeglab2022.1');

% Directory where your EEG datasets are stored
dataDir = 'C:/Users/juhoffmann/Desktop/EEG_BIDS/EEG_250Hz/Matlab';

% List of participants - update this list according to your participant IDs/names
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

participants = subjects(idx);      

% Initialize a structure to hold the aggregated N170 amplitudes and trial numbers
results = struct;

% Initialize variables for aggregated data
aggregatedN170Amplitudes = [];
aggregatedTrialNumbers = [];

for i = 1:length(participants)
    
    FileName = participants(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath); 
    participant = strcat('sub',FileName(5:7));   
        
    % Epoch extraction around the N170 component
    EEG = pop_epoch(EEG, {  'h_h_strong'  'h_h_weak'  'h_n_strong'  'h_n_weak'  'h_s_strong'  'h_s_weak'  'n_h_strong'  'n_h_weak'  'n_n_strong'  'n_n_weak'  'n_s_strong'  'n_s_weak'  's_h_strong'  's_h_weak'  's_n_strong'  's_n_weak'  's_s_strong'  's_s_weak'  },...
    [-0.2  0.8],  'epochinfo', 'yes');
    
    % Baseline correction
    EEG = pop_rmbase(EEG, [-200, 0]);
    
    % Specify the channels of interest and time window for N170
    chanIndices = find(ismember({EEG.chanlocs.labels}, {'P7', 'P8', 'PO7', 'PO8'}));
    timeWindow = [150, 200]; % N170 time window in milliseconds
    
    % Find the time indices for the N170 window
    timeIndices = find(EEG.times >= timeWindow(1) & EEG.times <= timeWindow(2));
    
    % Initialize variable to store N170 amplitudes for each trial
    n170Amplitudes = zeros(1, EEG.trials);
    
    % Loop through each trial to calculate N170 amplitude
    for t = 1:EEG.trials
        % Calculate mean amplitude across specified channels and time window
        n170Amplitudes(t) = mean(mean(EEG.data(chanIndices, timeIndices, t), 2), 1);
    end

    % Define trialNumbers here
    trialNumbers = 1:EEG.trials;

    % Append the participant's data to the aggregated variables
    aggregatedN170Amplitudes = [aggregatedN170Amplitudes, n170Amplitudes];
    aggregatedTrialNumbers = [aggregatedTrialNumbers, 1:EEG.trials]; % This assumes trial numbers restart at 1 for each participant
    
    % Store the N170 amplitudes for this participant
    results.(participant).n170_amplitudes = n170Amplitudes;
    
    % Analysis of N170 amplitude vs. trial number
    [b,~,~,~,stats] = regress(n170Amplitudes', [ones(EEG.trials, 1) trialNumbers']);
    
    % Store regression results for the individual participant
    results.(participant).regression.b = b;
    results.(participant).regression.stats = stats;
end

% Perform regression analysis on the aggregated data
[b_agg,~,~,~,stats_agg] = regress(aggregatedN170Amplitudes', [ones(length(aggregatedN170Amplitudes), 1), aggregatedTrialNumbers']);

% Store the aggregated regression results
results.aggregated.regression.b = b_agg;
results.aggregated.regression.stats = stats_agg;

%%
figure; % Create a new figure
scatter(aggregatedTrialNumbers, aggregatedN170Amplitudes, 'filled'); % Plot the aggregated data points
hold on; % Hold on to plot the regression line on the same figure

% Calculate the regression line
regLine = results.aggregated.regression.b(1) + results.aggregated.regression.b(2) * aggregatedTrialNumbers;

% Plot the regression line
plot(aggregatedTrialNumbers, regLine, 'r-', 'LineWidth', 2); % 'r-' for a red line
title('Aggregated N170 Amplitudes vs. Trial Numbers');
xlabel('Trial Number');
ylabel('N170 Amplitude (\muV)');
legend('Data points', 'Regression Line', 'Location', 'best');
hold off; % Release the figure for other plots


