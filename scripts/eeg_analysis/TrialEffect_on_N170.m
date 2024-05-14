% Add EEGLAB to MATLAB path (update this path to your EEGLAB installation)
addpath('C:/Users/juhoffmann/Desktop/eeglab2022.1');
eeglab;
% Directory where your EEG datasets are stored
DataPath = 'C:/Users/juhoffmann/Desktop/EEG_BIDS/EEG_250Hz/Matlab';

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
aggregatedN170Amplitudes =[];
aggregatedTrialNumbers =[];
% Define the conditions
conditions = {'h_h_strong' 'h_h_weak' 'h_n_strong' 'h_n_weak' 'h_s_strong' 'h_s_weak' ...
              'n_h_strong' 'n_h_weak' 'n_n_strong' 'n_n_weak' 'n_s_strong' 'n_s_weak' ...
              's_h_strong' 's_h_weak' 's_n_strong' 's_n_weak' 's_s_strong' 's_s_weak'};

% Loop through each participant
for i = 1:length(participants)
    
    FileName = participants(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath); 
    participant = strcat('sub', FileName(5:7));   
        
    % Epoch extraction around the N170 component for specific conditions only
    EEG = pop_epoch(EEG, conditions, [-0.2 0.8], 'epochinfo', 'yes');
    
    % Baseline correction
    EEG = pop_rmbase(EEG, [-200 0]);
    
    % Specify the channels of interest and time window for N170
    chanIndices = find(ismember({EEG.chanlocs.labels}, {'P7', 'P8', 'PO7', 'PO8'}));
    timeWindow = [150, 200]; % N170 time window in milliseconds
    timeIndices = find(EEG.times >= timeWindow(1) & EEG.times <= timeWindow(2));
    
    % Initialize variable to store N170 amplitudes for each trial
    n170Amplitudes = [];
    for j = 1:length(conditions)
        condition = conditions{j};
            % Find Epochs
            Epoch = extractfield(EEG.event,'epoch');

            % Find stimulus Type
            Type = extractfield(EEG.event,'type');

            % Find trials for the current condition
            conditionTrials = find(strcmp(Type, condition));
    
            conditionEpochs = Epoch(conditionTrials);

    % Access and calculate N170 amplitude only for specific conditions


        % Calculate mean amplitude for each of these trials
        for t = conditionEpochs
            mean_erp_range = mean(mean(EEG.data(chanIndices,timeIndices,t),3),1);
            n170Amplitude = min(mean_erp_range);
            n170Amplitudes = [n170Amplitudes, n170Amplitude];
        end
    end

    % Append the participant's data to the aggregated variables
    aggregatedN170Amplitudes = [aggregatedN170Amplitudes, n170Amplitudes];
    aggregatedTrialNumbers = [aggregatedTrialNumbers, 1:length(n170Amplitudes)]; % This assumes trial numbers restart at 1 for each participant
    
    % Store the N170 amplitudes and regression analysis results for this participant
    results.(participant).n170_amplitudes = n170Amplitudes;
    [b, ~, ~, ~, stats] = regress(n170Amplitudes', [ones(length(n170Amplitudes), 1) (1:length(n170Amplitudes))']);
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


%% get individual slope

% Initialize arrays to store participant IDs and their corresponding slopes
participantIDs = strings(length(participants), 1);
slopes = zeros(length(participants), 1);

for i = 1:length(participants)
    participant = strcat('sub', participants(i).name(5:7)); % Extract participant ID
    participantIDs(i) = participant; % Store participant ID
    slopes(i) = results.(participant).regression.b(2); % Store the slope of N170 amplitude change
end

% Create a table with participant IDs and their slopes
slopeTable = table(participantIDs, slopes, 'VariableNames', {'ParticipantID', 'Slope'});

% Display the table
disp(slopeTable);

% Optionally, save the table to a file
writetable(slopeTable, 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/N170_ParticipantSlopes.csv');

%% combine individual slope with existing datasheet

% Step 1: Read the existing dataset
existingDatasetPath = 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/erp_N170_HC_and_MDD.xlsx';
existingDataset = readtable(existingDatasetPath);

% Step 2: Read the slopes data
slopeDatasetPath = 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/N170_ParticipantSlopes.csv';
slopeDataset = readtable(slopeDatasetPath);

% Adjust participant IDs in the slope dataset to match the existing dataset format if necessary
% For example, if existing dataset uses "sub-01" format and slope uses "sub01"
slopeDataset.ParticipantID = regexprep(slopeDataset.ParticipantID, 'sub', 'sub-');

% Step 3: Match participant IDs and add slopes to the existing dataset
% Initialize a column for slopes in the existing dataset
existingDataset.Slope = NaN(height(existingDataset), 1); % Add a new 'Slope' column

for i = 1:height(slopeDataset)
    participantID = slopeDataset.ParticipantID{i};
    slope = slopeDataset.Slope(i);
    
    % Find the row index(es) for this participant in the existing dataset
    idx = find(strcmp(existingDataset.subName, participantID));
    
    % Step 4: Update the existing dataset with the slope information
    if ~isempty(idx)
        existingDataset.Slope(idx) = slope;
    end
end

% Step 5: Save the updated dataset back to an Excel file
updatedDatasetPath = 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/updated_erp_N170_HC_and_MDD.xlsx';
writetable(existingDataset, updatedDatasetPath);
