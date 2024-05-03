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


% Loop through each participant
for i = 1:length(participants)
    FileName = participants(i).name;
    filePath = fullfile(DataPath, FileName);
    EEG = pop_loadbva(filePath); 
    participant = strcat('sub',FileName(5:7)); 

    % Epoch extraction for all conditions at once
    conditions = {'h_h_strong' 'h_h_weak' 'h_n_strong' 'h_n_weak' 'h_s_strong' 'h_s_weak' 'n_h_strong' 'n_h_weak' 'n_n_strong' 'n_n_weak' 'n_s_strong' 'n_s_weak' 's_h_strong' 's_h_weak' 's_n_strong' 's_n_weak' 's_s_strong' 's_s_weak'};
    EEG = pop_epoch(EEG, conditions, [-0.2  0.8], 'epochinfo', 'yes');
    EEG = pop_rmbase(EEG, [-200, 0]);

    % Define the channels and time window for N170
    chanIndices = find(ismember({EEG.chanlocs.labels}, {'P7', 'P8', 'PO7', 'PO8'}));
    timeWindow = [150, 200]; % N170 time window in milliseconds
    timeIndices = find(EEG.times >= timeWindow(1) & EEG.times <= timeWindow(2));

    % Check indices are correct
    disp(['Channel Indices: ' num2str(chanIndices)]);
    disp(['Time Indices: ' num2str(timeIndices)]);

    % Initialize data structure for emotional categories
    emotionData = struct('Happy', [], 'Neutral', [], 'Sad', []);

    % Aggregate trials based on emotional category
    for j = 1:length(conditions)
        condition = conditions{j};
        emotionKey = ''; % Determine the emotion category
        if startsWith(condition, 'h')
            emotionKey = 'Happy';
        elseif startsWith(condition, 'n')
            emotionKey = 'Neutral';
        elseif startsWith(condition, 's')
            emotionKey = 'Sad';
        end

        % Find Epochs
        Epoch = extractfield(EEG.event,'epoch');

        % Find stimulus Type
        Type = extractfield(EEG.event,'type');

        % Find trials for the current condition
        conditionTrials = find(strcmp(Type, condition));
        
        conditionEpochs = Epoch(conditionTrials);

        % Debug: Check trial indices
        disp(['Condition: ' condition ' - Trials: ' num2str(conditionEpochs)])


        % Compute N170 amplitudes for these trials and aggregate them into the correct category
        for k = 1:length(conditionEpochs)
            t = conditionEpochs(k);
            if isempty(t) || isempty(chanIndices) || isempty(timeIndices)
                warning('Skipping empty indices for trial computation.');
                continue;
            end

            mean_erp_range = mean(mean(EEG.data(chanIndices,timeIndices,t),3),1);
            n170Amplitude = min(mean_erp_range);
            emotionData.(emotionKey) = [emotionData.(emotionKey), n170Amplitude];
        end
    end
    N170_values.(participant) = struct();
    N170_values.(participant) = emotionData;

    % Compute and store regression slopes for each emotional category
    results.(participant) = struct();
    results_trial_sets.(participant) = struct();
    emotions = fieldnames(emotionData);
    
    % Compute slopes for 10 trials per condition
    for emotion = emotions'
        n170Amplitudes = emotionData.(emotion{1});
        trialNumbers = 1:length(n170Amplitudes);
        trial_sets = floor(length(trialNumbers)/10); % Number of full 10-trial sets
        for k = 1:trial_sets
            % Select the trials for this set
            trial_indices = (k-1)*10 + (1:10);
            selected_trials = n170Amplitudes(trial_indices);
            % Compute slope
            [b, ~, ~, ~, stats] =  regress(selected_trials', [ones(length(selected_trials), 1) (1:length(selected_trials))']);
            results_trial_sets.(participant).(emotion{1})(k).slope = b(2);
            results_trial_sets.(participant).(emotion{1})(k).stats = stats;
        end
    
     % Compute slopes for every condition
     for emotion = emotions'
            n170Amplitudes = emotionData.(emotion{1});
            trialNumbers = 1:length(n170Amplitudes); % Assuming consecutive numbering, adjust if needed
            if isempty(n170Amplitudes)
                warning(['Skipping regression for empty N170 amplitudes in ' emotion{1}]);
                continue;
            end
            [b, ~, ~, ~, stats] = regress(n170Amplitudes', [ones(length(trialNumbers), 1) trialNumbers']);
            results.(participant).(emotion{1}).slope = b(2);
            results.(participant).(emotion{1}).stats = stats;
     end
    end
end





%% Plot N170 for every participant per condition
participants = fieldnames(N170_values); 
numParticipants = length(participants);

% Define number of emotions
emotions = {'Happy', 'Neutral', 'Sad'};
numEmotions = length(emotions);

% Create figure
figure;
for i = 1:numParticipants
    participant = participants{i};
    % Subplot for each participant
    for j = 1:numEmotions
        emotion = emotions{j};
        subplot(numParticipants, numEmotions, (i - 1) * numEmotions + j);
        plot(N170_values.(participant).(emotion));
        title([participant ' - ' emotion]);
        xlabel('Trial');
        ylabel('N170 Amplitude (µV)');
        grid on;
    end
end

% Adjust layout
sgtitle('N170 Values by Participant and Emotion');



%% Plot mean N170 for every condition
% Initialize storage for data aggregated across all participants
aggregateData = struct('Happy', [], 'Neutral', [], 'Sad', []);

participants = fieldnames(N170_values);
maxTrials = 0; % Find the maximum number of trials

% First determine the maximum number of trials to initialize vectors correctly
for i = 1:length(participants)
    partData = N170_values.(participants{i});
    for emotion = fieldnames(partData)'
        numTrials = length(partData.(emotion{1}));
        maxTrials = max(maxTrials, numTrials);
    end
end

% Initialize vectors for each emotion to maxTrials length
for emotion = fieldnames(aggregateData)'
    aggregateData.(emotion{1}) = cell(1, maxTrials);
end

% Loop through each participant's data to aggregate N170 values
for i = 1:length(participants)
    partData = N170_values.(participants{i});
    for emotion = fieldnames(partData)'
        emotionName = emotion{1};
        for t = 1:length(partData.(emotionName))
            aggregateData.(emotionName){t} = [aggregateData.(emotionName){t}, partData.(emotionName)(t)];
        end
    end
end

% Calculate the mean of each trial across all participants
meanData = struct('Happy', zeros(1, maxTrials), 'Neutral', zeros(1, maxTrials), 'Sad', zeros(1, maxTrials));
for emotion = fieldnames(meanData)'
    for t = 1:maxTrials
        % Only calculate the mean if there is data
        if ~isempty(aggregateData.(emotion{1}){t})
            meanData.(emotion{1})(t) = mean(aggregateData.(emotion{1}){t});
        end
    end
end


% Plotting
figure;
hold on;
colors = lines(3); % Get 3 distinct colors

% Plot each emotion's mean data
plot(meanData.Happy, 'Color', colors(1, :), 'DisplayName', 'Happy');
plot(meanData.Neutral, 'Color', colors(2, :), 'DisplayName', 'Neutral');
plot(meanData.Sad, 'Color', colors(3, :), 'DisplayName', 'Sad');

% Add labels and legend
xlabel('Trial Number');
ylabel('Mean N170 Amplitude (µV)');
title('Mean N170 Amplitude Across Trials by Emotion');
legend show;
grid on;



%% Extract Slopes

participantData = table();

% Collect data
for participant = fieldnames(results)'
    p = participant{1};
    % Check if all slopes are available
    if isfield(results.(p), 'Happy') && isfield(results.(p), 'Neutral') && isfield(results.(p), 'Sad')
        slope_happy = results.(p).Happy.slope;
        slope_neutral = results.(p).Neutral.slope;
        slope_sad = results.(p).Sad.slope;
        
        % Create a temporary table to append
        tempTable = table({p}, slope_happy, slope_neutral, slope_sad, ...
            'VariableNames', {'Participant', 'Slope_Happy', 'Slope_Neutral', 'Slope_Sad'});
        
        % Append to the main table
        participantData = [participantData; tempTable];
    else
        warning(['Missing data for participant ' p]);
    end
end

% Calculate mean slopes for each condition
meanSlopeHappy = mean(participantData.Slope_Happy);
meanSlopeNeutral = mean(participantData.Slope_Neutral);
meanSlopeSad = mean(participantData.Slope_Sad);

% Display the mean slopes
disp(['Mean Slope for Happy: ', num2str(meanSlopeHappy)]);
disp(['Mean Slope for Neutral: ', num2str(meanSlopeNeutral)]);
disp(['Mean Slope for Sad: ', num2str(meanSlopeSad)]);

allSlopes = [participantData.Slope_Happy; participantData.Slope_Neutral; participantData.Slope_Sad];

meanSlope = mean(allSlopes);

% Export to CSV file
writetable(participantData, 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/ParticipantSlopes_N170_emotions.csv');

%% Extract data in long format

participantDataLong = table();


% Collect data in long format with repetition
for participant = fieldnames(results)'
    p = participant{1};
    % Check if all slopes are available
    if isfield(results.(p), 'Happy') && isfield(results.(p), 'Neutral') && isfield(results.(p), 'Sad')
        slope_happy = results.(p).Happy.slope;
        slope_neutral = results.(p).Neutral.slope;
        slope_sad = results.(p).Sad.slope;
        
        % Create a temporary table to append (long format, repeated twice)
        emotions = {'Happy'; 'Sad'; 'Neutral'};
        slopes = [slope_happy; slope_sad; slope_neutral];
        tempTable = table(repmat({p}, 6, 1), ...
            repmat(slopes, 2, 1), ...
            repmat(emotions, 2, 1), ...
            'VariableNames', {'Participant', 'Slope', 'Emotion'});
        
        % Append to the main table
        participantDataLong = [participantDataLong; tempTable];
    else
        warning(['Missing data for participant ' p]);
    end
end


% Export to CSV file
writetable(participantDataLong, 'W:/Fmri_Forschung/Allerlei/JuliaS/GitHub/UnconsciousBias/data/ParticipantSlopesEmotionsLong.csv');


% 
% 
% 
% % Assuming slopes are in separate arrays: slopesHappy, slopesSad, slopesNeutral
% meanSlopeHappy = mean(slope_happy);
% meanSlopeSad = mean(slope_sad);
% meanSlopeNeutral = mean(slope_neutral);
% 
% % Display the mean slopes
% disp(['Mean Slope Happy: ', num2str(meanSlopeHappy)]);
% disp(['Mean Slope Sad: ', num2str(meanSlopeSad)]);
% disp(['Mean Slope Neutral: ', num2str(meanSlopeNeutral)]);

%%

% Define the emotions
emotions = {'Happy', 'Neutral', 'Sad'};

% Initialize variables to store slopes for plotting
plotData = struct('Happy', [], 'Neutral', [], 'Sad', []);

% Gather data from results_trial_sets
participants = fieldnames(results_trial_sets);
for i = 1:length(participants)
    participant = participants{i};
    for emotion = emotions
        emotionName = emotion{1};
        if isfield(results_trial_sets.(participant), emotionName)
            for k = 1:length(results_trial_sets.(participant).(emotionName))
                if isfield(results_trial_sets.(participant).(emotionName)(k), 'slope')
                    plotData.(emotionName) = [plotData.(emotionName), results_trial_sets.(participant).(emotionName)(k).slope];
                end
            end
        end
    end
end

% Plotting
figure;
hold on;

% Colors for each emotion
colors = lines(numel(emotions));  % Get distinct colors for each emotion

% Plot each emotion's data
for i = 1:numel(emotions)
    emotionName = emotions{i};
    plot(plotData.(emotionName), '-o', 'Color', colors(i, :), 'DisplayName', emotionName);
end

% Enhance the plot with labels, legend, and grid
xlabel('Set Number (each set represents 10 trials)');
ylabel('Slope');
title('Slopes of N170 Amplitudes Across Trial Sets for Each Emotion');
legend show;
grid on;

% Check if any data was empty and print a message if so
for emotion = emotions
    if isempty(plotData.(emotion{1}))
        fprintf('No data available for %s\n', emotion{1});
    end
end


%% plot per emotion

% Colors for each emotion, defined outside the loop for consistency
colors = lines(numel(emotions));  % Get distinct colors for each emotion

% Plot each emotion's data in separate figures
for i = 1:numel(emotions)
    emotionName = emotions{i};
    
    % Create a new figure for each emotion
    figure;
    hold on;
    
    % Plot the data
    plot(plotData.(emotionName), '-o', 'Color', colors(i, :), 'DisplayName', emotionName);
    
    % Enhance the plot with labels, legend, and grid
    xlabel('Set Number (each set represents 10 trials)');
    ylabel('Slope');
    title(['Slopes of N170 Amplitudes Across Trial Sets for ', emotionName]);
    legend show;
    grid on;
    
    % Check if any data was empty and print a message if so
    if isempty(plotData.(emotionName))
        fprintf('No data available for %s\n', emotionName);
    end
    
    hold off;
end

