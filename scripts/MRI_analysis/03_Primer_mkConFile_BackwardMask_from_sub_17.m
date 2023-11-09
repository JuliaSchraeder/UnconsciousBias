
function [names, onsets, durations] = Primer_mkConFile_BackwardMask_from_sub_17(subName)
root = '/bif/storage/storage1/projects/emocon/Data/';                       %Path to the data
projectfolder = '/bif/storage/storage1/projects/emocon';                    %Path to the projectfolder
logDir = fullfile(root,'Behav_data','Logfiles');                            %Path to the csv data

subjects = dir(fullfile(logDir,'*.csv'));                                   %find subjects with csv file generated during the PsychoPy experiment
subs = length(subjects);                                                    %get number of subjects
csv_names = {subjects(1:subs).name};                                        %extract csv file names
csv_names = csv_names';                                                     %convert 1x3 to 3x1

subNumber = extractAfter(subName, 3);                                       %get number of subjects
Index = find(contains(csv_names,subNumber));                                %find row number of subName in "csv_names"

 
logName = csv_names(Index);                                                 %get specific name of csv file for subName
logName = char(logName);                                                    %convert cell to string


%%Output Directory festlegen  
if ~exist(fullfile(projectfolder, 'ConFileFolder_Primer',subName))                 %generate a outputfolder for the conditions file for every participant
    mkdir(fullfile(projectfolder, 'ConFileFolder_Primer',subName));
end
outputFolder = fullfile(projectfolder,'ConFileFolder_Primer',subName);             %define it as output folder

%% *** Dateneinlesen und codieren *** %%
%read in the csv files
fileName = fullfile(logDir,logName);
fileID = fopen(fileName,'r','n','UTF-8');
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s';

fseek(fileID, 3, 'bof'); %open the file
data = textscan(fileID, formatSpec, 'Delimiter', ';', 'TextType', 'string', 'HeaderLines' ,1, 'ReturnOnError', false, 'EndOfLine', '\r\n'); % delimiter for 61 ",", for other ";"
fclose(fileID);

%Get Onsets
MR_Onset_row = double(data{1,35});                                          % time in seconds for the first MR Trigger is in column 37 row 1
ISI_row = double(data{1,41});                                               % first fixation cross row
ISI_start = ISI_row(3,1);                                                   % the first trial starts in row 3 because the first two rows are example trials
MR_Onset = MR_Onset_row(~isnan(MR_Onset_row));                              % delete NaNs in this row
Primer_Onset_row = str2double(data{1,45});                                  % get times in seconds of target onset
Target_Onset_row = str2double(data{1,49});                                  % get times in seconds of target onset
PrimerTime_row = str2double(data{1,17});                                    % get duration of primer
PrimerTime_row = PrimerTime_row*(1/0.12);

%Get Stimulus Condition
strPrimer_Emotion = data{1,4};                                              % get primer emotion
strTarget_Emotion = data{1,9};                                              % get target emotion
strMask_Condition = data{1,17};                                             % get mask condition (strongly masked = 0.016 primer duration, weakly masked = 150 primer duration; mask is always 66.7ms long)

%Transform Stimulus Condition from string to number
Primer_Emotion = strrep(strPrimer_Emotion, 'neutral', '1');                 %rename primer neutral = 1
Primer_Emotion = strrep(Primer_Emotion, 'sad', '2');                        %rename primer sad = 2
Primer_Emotion = strrep(Primer_Emotion, 'happy', '3');                      %rename primer happy = 3
Primer_Emotion = str2double(Primer_Emotion);                                %transform Primer information sting to double

Target_Emotion = strrep(strTarget_Emotion, 'neutral', '1');                 %rename target neutral = 1
Target_Emotion = strrep(Target_Emotion, 'sad', '2');                        %rename target sad = 2
Target_Emotion = strrep(Target_Emotion, 'happy', '3');                      %rename target happy = 3
Target_Emotion = str2double(Target_Emotion);                                %transform Target information sting to double

Mask_Condition = str2double(strMask_Condition);                             %transform Mask condition information string to double,           


%Delete first two example rows and NaNs in last row:
Primer_Onset_row([1,2,363],:) = [];
Target_Onset_row([1,2,363],:) = [];
PrimerTime_row([1,2,363],:) = [];
Mask_Condition([1,2,363],:) = [];
Primer_Emotion([1,2,363],:) = [];
Target_Emotion([1,2,363],:) = [];

%Substract Onsets with MR Onset
Target_Onset_row = Target_Onset_row - MR_Onset -5.4;
Primer_Onset_row = Primer_Onset_row - MR_Onset -5.4;


%% Get Onsets %%

%Condition onset is onset of target. Onsets are taken from the
%Target_Onset_row and depend on the primer, target and mask condition
%Mask Condition = 2 --> strongly masked trial (0.016ms = 2 frames)
%Mask Condition = 18 --> weakly masked trial (150ms = 18 frames)
%prime emotion = 3 and target emotion = 3 and mask contion = 2 --> happy_happy_unconscious trials

onsets_happy_unconscious = Primer_Onset_row(Primer_Emotion == 3 & Mask_Condition == 2);
onsets_sad_unconscious = Primer_Onset_row(Primer_Emotion == 2 & Mask_Condition == 2);
onsets_neutral_unconscious = Primer_Onset_row(Primer_Emotion == 1 & Mask_Condition == 2);

onsets_happy_conscious = Primer_Onset_row(Primer_Emotion == 3 & Mask_Condition == 18);
onsets_sad_conscious = Primer_Onset_row(Primer_Emotion == 2 & Mask_Condition == 18);
onsets_neutral_conscious = Primer_Onset_row(Primer_Emotion == 1 & Mask_Condition == 18);

%onsets_primer = Primer_Onset_row;
onsets_intro = 0;

%% Get Durations %%

duration_happy_unconscious = zeros(length(onsets_happy_unconscious),1);     %Duration is zero
duration_sad_unconscious = zeros(length(onsets_sad_unconscious),1); 
duration_neutral_unconscious = zeros(length(onsets_neutral_unconscious),1); 

duration_happy_conscious = zeros(length(onsets_happy_conscious),1); 
duration_sad_conscious = zeros(length(onsets_sad_conscious),1); 
duration_neutral_conscious = zeros(length(onsets_neutral_conscious),1); 

%duration_primer = PrimerTime_row/1000;                                     %PrimerTime_row is in ms, not in seconds. Therefore divide it with 1000
duration_intro = ISI_start - MR_Onset -5.4;                                      %substract the start of the experiment with MR onset to get the duration of the instruction


%% Create cells for Conditionsfile

names{1} = 'happy_unconscious';                                      
names{2} = 'sad_unconscious';
names{3} = 'neutral_unconscious';
names{4} = 'happy_conscious';
names{5} = 'sad_conscious';
names{6} = 'neutral_conscious'; 
names{7} = 'intro';                                                       


onsets{1} = onsets_happy_unconscious;                                 
onsets{2} = onsets_sad_unconscious;
onsets{3} = onsets_neutral_unconscious;
onsets{4} = onsets_happy_conscious;
onsets{5} = onsets_sad_conscious;
onsets{6} = onsets_neutral_conscious;
onsets{7} = onsets_intro;                                                 


durations{1} = duration_happy_unconscious;                            
durations{2} = duration_sad_unconscious;
durations{3} = duration_neutral_unconscious;
durations{4} = duration_happy_conscious;
durations{5} = duration_sad_conscious;
durations{6} = duration_neutral_conscious;
durations{7} = duration_intro;    

save(fullfile(outputFolder,'ConFile_BackwardMask'), 'names', 'onsets','durations');
end

