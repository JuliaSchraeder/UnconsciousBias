clear

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

% Initialize data structures

conditions = {'h_h_strong', 'h_h_weak', 'h_n_strong', 'h_n_weak', 'h_s_strong', 'h_s_weak', ...
              'n_h_strong', 'n_h_weak', 'n_n_strong', 'n_n_weak', 'n_s_strong', 'n_s_weak', ...
              's_h_strong', 's_h_weak', 's_n_strong', 's_n_weak', 's_s_strong', 's_s_weak'};
emotionCategories = {'h', 'n', 's'};
emotionLabels = {'Happy', 'Neutral', 'Sad'};
setSize = 10; % Trials per set

% Initialize results structure
results = struct;

% Loop through each participant
for i = 1:length(participants)
    FileName = participants(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath);
    participant = strcat('sub', FileName(5:7));
    results.(participant) = struct;

    % Epoch extraction and baseline correction
    EEG = pop_epoch(EEG, conditions, [-0.2 0.8], 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-200 0]);

    % Channels and time window for N170
    chanIndices = find(ismember({EEG.chanlocs.labels}, {'P7', 'P8', 'PO7', 'PO8'}));
    timeWindow = [150, 200];
    timeIndices = find(EEG.times >= timeWindow(1) & EEG.times <= timeWindow(2));

    % Process each condition
    for emotion = 1:length(emotionLabels)
        emotionKey = emotionLabels{emotion};
        emotionData = [];
        
        % Collect data for the emotion
        for condition = conditions(startsWith(conditions, emotionCategories{emotion}))
            Epoch = extractfield(EEG.event, 'epoch');
            Type = extractfield(EEG.event, 'type');
            conditionTrials = find(strcmp(Type, condition));
            conditionEpochs = Epoch(conditionTrials);

            for t = conditionEpochs
                minAmplitude = min(mean(EEG.data(chanIndices, timeIndices, t), 1));
                emotionData = [emotionData, minAmplitude];
            end
        end
        
        % Store emotion-specific data
        results.(participant).(emotionKey).data = emotionData;

        % Calculate slope for all trials
        if ~isempty(emotionData)
            [b, ~, ~, ~, stats] = regress(emotionData', [ones(length(emotionData), 1) (1:length(emotionData))']);
            results.(participant).(emotionKey).slope = b(2);  % Slope of regression
            results.(participant).(emotionKey).stats = stats; % Regression stats
        end

        % Calculate slope for every set of 10 trials
        numSets = floor(length(emotionData) / setSize);
        setSlopes = [];
        for setNum = 1:numSets
            startIndex = (setNum-1) * setSize + 1;
            endIndex = startIndex + setSize - 1;
            setAmplitudes = emotionData(startIndex:endIndex);
            [b, ~, ~, ~, ~] = regress(setAmplitudes', [ones(size(setAmplitudes')), (1:length(setAmplitudes))']);
            setSlopes = [setSlopes, b(2)];
        end
        results.(participant).(emotionKey).setSlopes = setSlopes;
    end
end


%% Plot
clear
% Load the data
load('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/scripts/eeg_analysis/N170results.mat');

% Initialize variables for plotting
emotions = {'Happy', 'Sad', 'Neutral'};
participantNames = fieldnames(results); % Get all participant names
numTrials = 120; % Assuming a maximum of 50 trials per participant for each emotion

% Prepare to store the mean values and slopes
meanN170 = struct();
regressionSlopes = zeros(1, numel(emotions)); % Store the slopes of regression lines

% Initialize data storage for each emotion and trial
for e = 1:numel(emotions)
    emotion = emotions{e};
    meanN170.(emotion) = zeros(1, numTrials); % Preallocate with zeros
    countTrials = zeros(1, numTrials); % To count the number of participants per trial
end

% Collect N170 values across all participants
for p = 1:numel(participantNames)
    participant = participantNames{p};
    for e = 1:numel(emotions)
        emotion = emotions{e};
        if isfield(results.(participant), emotion) && ~isempty(results.(participant).(emotion).data)
            dataLength = length(results.(participant).(emotion).data);
            meanN170.(emotion)(1:dataLength) = meanN170.(emotion)(1:dataLength) + results.(participant).(emotion).data;
            countTrials(1:dataLength) = countTrials(1:dataLength) + 1;
        end
    end
end

% Calculate the mean N170 values
for e = 1:numel(emotions)
    emotion = emotions{e};
    meanN170.(emotion)(countTrials > 0) = meanN170.(emotion)(countTrials > 0) ./ countTrials(countTrials > 0);
end

% Plot the results with regression lines
figure;
hold on;
colors = [0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330]; % Custom colors: less bright
legends = {};
for e = 1:numel(emotions)
    emotion = emotions{e};
    trials = 1:length(meanN170.(emotion)(meanN170.(emotion) ~= 0)); % Only consider trials with data
    meanValues = meanN170.(emotion)(meanN170.(emotion) ~= 0);
    plot(trials, meanValues, 'o', 'Color', colors(e,:), 'DisplayName', emotion); % Display as dots
    legends{end+1} = emotion; % Add to legend
    
    % Fit linear regression
    mdl = fitlm(trials', meanValues');
    coef = mdl.Coefficients.Estimate;
    xFit = [min(trials), max(trials)];
    yFit = coef(1) + coef(2) * xFit;
    plot(xFit, yFit, '-', 'Color', colors(e,:), 'LineWidth', 2); % Display as lines
    
    % Store the slope
    regressionSlopes(e) = coef(2);
end

% Enhance the plot
xlabel('Trial Number');
ylabel('Mean N170 Amplitude');
title('Mean N170 Amplitude Across All Participants for Each Emotion with Regression Lines');

% Display the slopes for each emotion
for e = 1:numel(emotions)
    disp(['Mean slope for ', emotions{e}, ': ', num2str(regressionSlopes(e))]);
end


legends = {'happy', 'happy slope', 'sad','sad slope','neutral','neutral slope'};
% Create custom legend with matching colors
lgd = legend(legends, 'Location', 'NorthEastOutside'); legend('boxoff') % Adjust the legend location
hold off;



%% Plot multiple regression lines (every 10 trials)

% Load the data
load('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/scripts/eeg_analysis/N170results.mat');

% Initialize variables for plotting
emotions = {'Happy', 'Sad', 'Neutral'};
participantNames = fieldnames(results); % Get all participant names
numTrials = 120; % Total number of trials

% Prepare to store the mean values and slopes
meanN170 = struct();
regressionSlopes = zeros(1, numel(emotions)); % Store the slopes of regression lines

% Initialize data storage for each emotion and trial
for e = 1:numel(emotions)
    emotion = emotions{e};
    meanN170.(emotion) = zeros(1, numTrials); % Preallocate with zeros
    countTrials = zeros(1, numTrials); % To count the number of participants per trial
end

% Collect N170 values across all participants
for p = 1:numel(participantNames)
    participant = participantNames{p};
    for e = 1:numel(emotions)
        emotion = emotions{e};
        if isfield(results.(participant), emotion) && ~isempty(results.(participant).(emotion).data)
            dataLength = length(results.(participant).(emotion).data);
            meanN170.(emotion)(1:dataLength) = meanN170.(emotion)(1:dataLength) + results.(participant).(emotion).data;
            countTrials(1:dataLength) = countTrials(1:dataLength) + 1;
        end
    end
end

% Calculate the mean N170 values
for e = 1:numel(emotions)
    emotion = emotions{e};
    meanN170.(emotion)(countTrials > 0) = meanN170.(emotion)(countTrials > 0) ./ countTrials(countTrials > 0);
end

% Plot the results with regression lines every 10 trials
figure;
hold on;
colors = [0.4660 0.6740 0.1880; 0.9290 0.6940 0.1250; 0.3010 0.7450 0.9330]; % Custom colors: less bright

for e = 1:numel(emotions)
    emotion = emotions{e};
    meanValues = meanN170.(emotion);
    trials = 1:length(meanValues);
    % Plot the mean N170 values (dots)
    plot(trials, meanValues, 'o', 'Color', colors(e,:), 'DisplayName', emotion); % Display as dots
    
    
    % Fit linear regression for every 10 trials
    for i = 1:10:length(trials)-9
        stack_trials = trials(i:i+9);
        stack_values = meanValues(i:i+9);
        % Fit linear regression for each stack of 10 trials
        mdl = fitlm(stack_trials', stack_values');
        coef = mdl.Coefficients.Estimate;
        xFit = [min(stack_trials), max(stack_trials)];
        yFit = coef(1) + coef(2) * xFit;
        plot(xFit, yFit, '-', 'Color', colors(e,:), 'LineWidth', 1); % Display as lines
        % Store the slope for each stack of 10 trials
        regressionSlopes(e,i) = coef(2);
    end
end

% Enhance the plot
xlabel('Trial Number');
ylabel('Mean N170 Amplitude');
title('Mean N170 Amplitude Across All Participants for Each Emotion with Multiple Regression Lines for Stacks of 10 Trials');

% Display the slopes for each emotion
for e = 1:numel(emotions)
    emotion = emotions{e};
    disp(['Mean slope for ', emotion, ': ', num2str(mean(regressionSlopes(e,:)))]);
end

% Create custom legend with matching colors
%legend(legends, 'Location', 'NorthEastOutside'); legend('boxoff'); % Adjust the legend location
hold off;

%% Extract data

% Load the data
load('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/scripts/eeg_analysis/N170results.mat');

% Initialize variables
emotions = {'Happy', 'Sad', 'Neutral'};
participantNames = fieldnames(results); % Get all participant names

% Prepare to store the data
participantData = {};
emotionData = {};
meanN170Data = [];
slopeData = [];

% Process each participant
for p = 1:numel(participantNames)
    participant = participantNames{p};
    % Process each emotion for the current participant
    for e = 1:numel(emotions)
        emotion = emotions{e};
        % Check if the current participant and emotion have data
        if isfield(results.(participant), emotion) && ~isempty(results.(participant).(emotion).data)
            data = results.(participant).(emotion).data;
            numTrials = size(data, 2); % Get the total number of trials
            % Determine the approximate number of stacks
            numStacks = ceil(numTrials / 10);
            % Determine the size of each stack
            stackSize = ceil(numTrials / numStacks);
            % Initialize arrays to store mean N170 and slope for each stack
            stackMeans = zeros(numStacks, 1);
            stackSlopes = zeros(numStacks, 1);
            % Iterate over each stack
            for s = 1:numStacks
                % Determine the start and end indices of the current stack
                startIndex = (s - 1) * stackSize + 1;
                endIndex = min(startIndex + stackSize - 1, numTrials);
                % Extract the N170 data for the current stack
                stackData = data(:, startIndex:endIndex);
                % Calculate the mean N170 for the current stack
                stackMeans(s) = mean(stackData, 'all');
                % Calculate the slope for the current stack
                trials = 1:numel(stackData);
                mdl = fitlm(trials', stackData');
                stackSlopes(s) = mdl.Coefficients.Estimate(2);
            end
            % Store the data in long format
            participantData = [participantData; repmat({participant}, numStacks, 1)];
            emotionData = [emotionData; repmat({emotion}, numStacks, 1)];
            meanN170Data = [meanN170Data; stackMeans(:)];
            slopeData = [slopeData; stackSlopes(:)];
        end
    end
end

% Create table
T = table(participantData, emotionData, meanN170Data, slopeData, 'VariableNames', {'Participant', 'Emotion', 'MeanN170', 'Slope'});

% Write to Excel
filename = 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/N170_Slope_Data_VariableStacks.xlsx';
writetable(T, filename);
disp(['Data exported to ' filename]);


%% Merge data with existing sheet

% Load the ERP dataset
erp_data = readtable('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/erp_N170_HC_and_MDD.xlsx');

% Load the slope dataset
slope_data = readtable('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/N170_Slope_Data_VariableStacks.xlsx');


% Normalize the 'subName' column in the ERP data
erp_data.subName = strrep(erp_data.subName, '-', '');

% Normalize the 'emotion' column in the ERP data
erp_data.emotion = cellfun(@(x) [upper(x(1)), lower(x(2:end))], erp_data.emotion, 'UniformOutput', false);

% Normalize the 'Participant' column in the slope data
slope_data.Participant = cellfun(@(x) strrep(x, '-', ''), slope_data.Participant, 'UniformOutput', false);

% Normalize the 'Emotion' column in the slope data
slope_data.Emotion = cellfun(@(x) [upper(x(1)), lower(x(2:end))], slope_data.Emotion, 'UniformOutput', false);

% Make the column names consistent
slope_data.Properties.VariableNames{'Participant'} = 'subName';
slope_data.Properties.VariableNames{'Emotion'} = 'emotion';

% Merge the datasets
merged_data = outerjoin(erp_data, slope_data, 'Keys', {'subName', 'emotion'}, 'MergeKeys', true, 'Type', 'inner');

% Save the merged dataset to a new Excel file
writetable(merged_data, 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/merged_N170_datasets_matlab.xlsx');










%% Test

% Display the structure of the 'results' matrix
disp(fieldnames(results));  % This will show the top-level fields (participant identifiers)

% Assuming 'participant1' is a valid field, adjust as necessary
disp(fieldnames(results.sub004));  % This will show the fields for a specific participant

% Display an example of stored data for a specific emotion under a participant
if isfield(results.sub004, 'happy')
    disp(results.sub004.happy)
end


% %%
% participants = subjects(idx);   
% 
% % Initialize data structures
% conditions = {'h_h_strong', 'h_h_weak', 'h_n_strong', 'h_n_weak', 'h_s_strong', 'h_s_weak', ...
%               'n_h_strong', 'n_h_weak', 'n_n_strong', 'n_n_weak', 'n_s_strong', 'n_s_weak', ...
%               's_h_strong', 's_h_weak', 's_n_strong', 's_n_weak', 's_s_strong', 's_s_weak'};
% emotionCategories = {'h', 'n', 's'};
% emotionLabels = {'Happy', 'Neutral', 'Sad'};
% setSize = 10; % Trials per set
% 
% % Initialize results structure
% results = struct;
% aggregatedN170Amplitudes = struct('Happy', [], 'Neutral', [], 'Sad', []);
% aggregatedTrialNumbers = struct('Happy', [], 'Neutral', [], 'Sad', []);
% 
% % Loop through each participant
% for i = 1:length(participants)
%     FileName = participants(i).name;
%     filePath = fullfile(DataPath, FileName);
%     EEG = pop_loadbva(filePath);
%     participant = strcat('sub', FileName(5:7));
%     results.(participant) = struct;
% 
%     % Epoch extraction and baseline correction
%     EEG = pop_epoch(EEG, conditions, [-0.2 0.8], 'epochinfo', 'yes');
%     EEG = pop_rmbase(EEG, [-200 0]);
% 
%     % Channels and time window for N170
%     chanIndices = find(ismember({EEG.chanlocs.labels}, {'P7', 'P8', 'PO7', 'PO8'}));
%     timeWindow = [150, 200];
%     timeIndices = find(EEG.times >= timeWindow(1) & EEG.times <= timeWindow(2));
% 
%     % Process each condition
%     for emotion = 1:length(emotionLabels)
%         emotionKey = emotionLabels{emotion};
%         emotionData = [];
%         
%         % Collect data for the emotion
%         for condition = conditions(startsWith(conditions, emotionCategories{emotion}))
%             Epoch = extractfield(EEG.event, 'epoch');
%             Type = extractfield(EEG.event, 'type');
%             conditionTrials = find(strcmp(Type, condition));
%             conditionEpochs = Epoch(conditionTrials);
% 
%             for t = conditionEpochs
%                 minAmplitude = min(mean(EEG.data(chanIndices, timeIndices, t), 1));
%                 emotionData = [emotionData, minAmplitude];
%             end
%         end
%         
%         % Store emotion-specific data
%         results.(participant).(emotionKey).data = emotionData;
% 
%         % Append the participant's data to the aggregated variables
%         aggregatedN170Amplitudes.(emotionKey) = [aggregatedN170Amplitudes.(emotionKey), emotionData];
%         aggregatedTrialNumbers.(emotionKey) = [aggregatedTrialNumbers.(emotionKey), 1:length(emotionData)];
%     end
% end
% 
% % Perform regression analysis on the aggregated data for each emotion
% for emotion = emotionLabels
%     key = emotion{1};
%     [b_agg,~,~,~,stats_agg] = regress(aggregatedN170Amplitudes.(key)', [ones(length(aggregatedN170Amplitudes.(key)), 1), aggregatedTrialNumbers.(key)']);
%     
%     % Store the aggregated regression results
%     results.aggregated.(key).regression.b = b_agg;
%     results.aggregated.(key).regression.stats = stats_agg;
% end
% 
% 
% 
% %%
% 
% load('W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/scripts/eeg_analysis/N170results_with_aggregrated_N170_per_emotion.mat');
% 
% 
% % Create a new figure
% figure;
% hold on; % Hold on to plot multiple datasets
% 
% light_colors = {[0.7, 0.7, 1], [0.7, 1, 0.7], [1, 0.7, 0.7]}; % Light blue, light green, light red
% dark_colors = {'b', 'g', 'r'}; % Darker colors for "Happy", "Neutral", "Sad"
% 
% y_min = inf;
% y_max = -inf;
% 
% % First pass: Find global y-axis limits
% for emotionIdx = 1:length(emotionLabels)
%     key = emotionLabels{emotionIdx};
%     aggregatedN170Amplitudes_emotion = aggregatedN170Amplitudes.(key);
%     y_min = min(y_min, min(aggregatedN170Amplitudes_emotion));
%     y_max = max(y_max, max(aggregatedN170Amplitudes_emotion));
% end
% 
% % Second pass: Plot with consistent y-axis limits
% for emotionIdx = 1:length(emotionLabels)
%     key = emotionLabels{emotionIdx};
%     light_color = light_colors{emotionIdx};
%     dark_color = dark_colors{emotionIdx};
% 
%     % Extract the data
%     aggregatedN170Amplitudes_emotion = aggregatedN170Amplitudes.(key);
%     aggregatedTrialNumbers_emotion = aggregatedTrialNumbers.(key);
%     b_agg = results.aggregated.(key).regression.b;
% 
%     % Create a subplot for each emotion
%     subplot(1, 3, emotionIdx);
%     hold on;
% 
%     % Plot the aggregated data points with some transparency
%     scatter(aggregatedTrialNumbers_emotion, aggregatedN170Amplitudes_emotion, 36, ...
%         'MarkerFaceColor', light_color, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
% 
%     % Calculate the regression line
%     regLine = b_agg(1) + b_agg(2) * aggregatedTrialNumbers_emotion;
% 
%     % Plot the regression line
%     plot(aggregatedTrialNumbers_emotion, regLine, [dark_color, '-'], 'LineWidth', 2);
% 
%     % Set plot titles and labels
%     title(key);
%     xlabel('Trial Number');
%     ylabel('N170 Amplitude (\muV)');
%     ylim([y_min y_max]); % Set consistent y-axis limits
%     legend('Data Points', 'Regression Line', 'Location', 'best');
% 
%     hold off;
% end
% 
% % Adjust layout
% sgtitle('Aggregated N170 Amplitudes vs. Trial Numbers for Different Emotions');
% 
% %% plot all 4 regressions
% 
% figure;
% 
% % Colors
% light_colors = {[0.7, 0.7, 1], [0.7, 1, 0.7], [1, 0.7, 0.7]};
% dark_colors = {'b', 'g', 'r'};
% 
% % First, we will find the common y-axis limits for all plots
% y_min = inf;
% y_max = -inf;
% allAmplitudes = [];
% allTrialNumbers = [];
% 
% for emotionIdx = 1:length(emotionLabels)
%     key = emotionLabels{emotionIdx};
%     aggregatedN170Amplitudes_emotion = aggregatedN170Amplitudes.(key);
%     y_min = min(y_min, min(aggregatedN170Amplitudes_emotion));
%     y_max = max(y_max, max(aggregatedN170Amplitudes_emotion));
%     allAmplitudes = [allAmplitudes, aggregatedN170Amplitudes_emotion];
%     allTrialNumbers = [allTrialNumbers, aggregatedTrialNumbers.(key)];
% end
% 
% % Include the combined amplitudes in the y-axis calculation
% y_min = min(y_min, min(allAmplitudes));
% y_max = max(y_max, max(allAmplitudes));
% 
% % Now, plot each emotion separately
% for emotionIdx = 1:length(emotionLabels)
%     key = emotionLabels{emotionIdx};
%     light_color = light_colors{emotionIdx};
%     dark_color = dark_colors{emotionIdx};
% 
%     % Extract the data
%     aggregatedN170Amplitudes_emotion = aggregatedN170Amplitudes.(key);
%     aggregatedTrialNumbers_emotion = aggregatedTrialNumbers.(key);
%     b_agg = results.aggregated.(key).regression.b;
% 
%     % Create a subplot for each emotion
%     subplot(2, 2, emotionIdx);
%     hold on;
% 
%     % Plot the aggregated data points with some transparency
%     scatter(aggregatedTrialNumbers_emotion, aggregatedN170Amplitudes_emotion, 36, ...
%         'MarkerFaceColor', light_color, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
% 
%     % Calculate the regression line
%     regLine = b_agg(1) + b_agg(2) * aggregatedTrialNumbers_emotion;
% 
%     % Plot the regression line
%     plot(aggregatedTrialNumbers_emotion, regLine, [dark_color, '-'], 'LineWidth', 2);
% 
%     % Set plot titles and labels
%     title(['N170 Amplitudes vs. Trial Numbers - ', key]);
%     xlabel('Trial Number');
%     ylabel('N170 Amplitude (\muV)');
%     ylim([y_min y_max]);
%     legend('Data Points', 'Regression Line', 'Location', 'best');
% 
%     hold off;
% end
% 
% % Combined plot with all data
% subplot(2, 2, 4);
% hold on;
% color = [0.5, 0.5, 0.5]; % Gray color for combined data
% 
% % Plot the aggregated data points with transparency
% scatter(allTrialNumbers, allAmplitudes, 36, ...
%     'MarkerFaceColor', color, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5);
% 
% % Calculate the combined regression line
% b_agg_all = results.aggregated.regression.b;
% regLine_all = b_agg_all(1) + b_agg_all(2) * allTrialNumbers;
% 
% % Plot the regression line
% plot(allTrialNumbers, regLine_all, 'k-', 'LineWidth', 2); % Black for the combined regression line
% 
% % Set plot titles and labels
% title('N170 Amplitudes vs. Trial Numbers - All Emotions');
% xlabel('Trial Number');
% ylabel('N170 Amplitude (\muV)');
% ylim([y_min y_max]);
% legend('Data Points', 'Regression Line', 'Location', 'best');
% 
% hold off;
% 
% % Adjust layout
% sgtitle('N170 Amplitudes vs. Trial Numbers');
% 
% 
% 
