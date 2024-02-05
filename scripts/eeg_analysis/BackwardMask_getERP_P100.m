function Erp_P100 = BackwardMask_getERP_P100(EEG,channel1,channel2,channel3,channel4,range_min,range_max)
  
 
% Create Epochs from -200ms to 800 ms
EEG = pop_epoch( EEG, {  'h_h_strong'  'h_h_weak'  'h_n_strong'  'h_n_weak'  'h_s_strong'  'h_s_weak'  'n_h_strong'  'n_h_weak'  'n_n_strong'  'n_n_weak'  'n_s_strong'  'n_s_weak'  's_h_strong'  's_h_weak'  's_n_strong'  's_n_weak'  's_s_strong'  's_s_weak'  },...
    [-0.2  0.8],  'epochinfo', 'yes');

% Remove baseline
EEG = pop_rmbase( EEG, [-200 0] ,[]);

% Find Epochs
Epoch = extractfield(EEG.event,'epoch');

% Find stimulus Type
Type = extractfield(EEG.event,'type');

%% Get Index of stimulus

% Bewusst
index_h_h_weak = find(strcmp(Type, 'h_h_weak'));
index_h_n_weak = find(strcmp(Type, 'h_n_weak'));
index_h_s_weak = find(strcmp(Type, 'h_s_weak'));

index_s_h_weak = find(strcmp(Type, 's_h_weak'));
index_s_n_weak = find(strcmp(Type, 's_n_weak'));
index_s_s_weak = find(strcmp(Type, 's_s_weak'));

index_n_h_weak = find(strcmp(Type, 'n_h_weak'));
index_n_n_weak = find(strcmp(Type, 'n_n_weak'));
index_n_s_weak = find(strcmp(Type, 'n_s_weak'));

% Unbewusst
index_h_h_strong = find(strcmp(Type, 'h_h_strong'));
index_h_n_strong = find(strcmp(Type, 'h_n_strong'));
index_h_s_strong = find(strcmp(Type, 'h_s_strong'));

index_s_h_strong = find(strcmp(Type, 's_h_strong'));
index_s_n_strong = find(strcmp(Type, 's_n_strong'));
index_s_s_strong = find(strcmp(Type, 's_s_strong'));

index_n_h_strong = find(strcmp(Type, 'n_h_strong'));
index_n_n_strong = find(strcmp(Type, 'n_n_strong'));
index_n_s_strong = find(strcmp(Type, 'n_s_strong'));


%% Get Trials with respective index
% Bewusst
trials_h_h_weak = Epoch(index_h_h_weak);
trials_h_n_weak = Epoch(index_h_n_weak);
trials_h_s_weak = Epoch(index_h_s_weak);

trials_s_h_weak = Epoch(index_s_h_weak);
trials_s_n_weak = Epoch(index_s_n_weak);
trials_s_s_weak = Epoch(index_s_s_weak);

trials_n_h_weak = Epoch(index_n_h_weak);
trials_n_n_weak = Epoch(index_n_n_weak);
trials_n_s_weak = Epoch(index_n_s_weak);

% Unbewusst
trials_h_h_strong = Epoch(index_h_h_strong);
trials_h_n_strong = Epoch(index_h_n_strong);
trials_h_s_strong = Epoch(index_h_s_strong);

trials_s_h_strong = Epoch(index_s_h_strong);
trials_s_n_strong = Epoch(index_s_n_strong);
trials_s_s_strong = Epoch(index_s_s_strong);

trials_n_h_strong = Epoch(index_n_h_strong);
trials_n_n_strong = Epoch(index_n_n_strong);
trials_n_s_strong = Epoch(index_n_s_strong);

%% Create a mean ERP for all trials 
% Bewusst
erp_epoch_h_h_weak = mean(EEG.data(:,:,trials_h_h_weak),3);
erp_epoch_h_n_weak = mean(EEG.data(:,:,trials_h_n_weak),3);
erp_epoch_h_s_weak = mean(EEG.data(:,:,trials_h_s_weak),3);

erp_epoch_s_h_weak = mean(EEG.data(:,:,trials_s_h_weak),3);
erp_epoch_s_n_weak = mean(EEG.data(:,:,trials_s_n_weak),3);
erp_epoch_s_s_weak = mean(EEG.data(:,:,trials_s_s_weak),3);

erp_epoch_n_h_weak = mean(EEG.data(:,:,trials_n_h_weak),3);
erp_epoch_n_n_weak = mean(EEG.data(:,:,trials_n_n_weak),3);
erp_epoch_n_s_weak = mean(EEG.data(:,:,trials_n_s_weak),3);

% Unbewusst
erp_epoch_h_h_strong = mean(EEG.data(:,:,trials_h_h_strong),3);
erp_epoch_h_n_strong = mean(EEG.data(:,:,trials_h_n_strong),3);
erp_epoch_h_s_strong = mean(EEG.data(:,:,trials_h_s_strong),3);

erp_epoch_s_h_strong = mean(EEG.data(:,:,trials_s_h_strong),3);
erp_epoch_s_n_strong = mean(EEG.data(:,:,trials_s_n_strong),3);
erp_epoch_s_s_strong = mean(EEG.data(:,:,trials_s_s_strong),3);

erp_epoch_n_h_strong = mean(EEG.data(:,:,trials_n_h_strong),3);
erp_epoch_n_n_strong = mean(EEG.data(:,:,trials_n_n_strong),3);
erp_epoch_n_s_strong = mean(EEG.data(:,:,trials_n_s_strong),3);

%% Create a mean for the Channels we are interested in
% Bewusst
mean_h_h_weak = mean(erp_epoch_h_h_weak([channel1,channel2,channel3,channel4],:),1);
mean_h_n_weak = mean(erp_epoch_h_n_weak([channel1,channel2,channel3,channel4],:),1);
mean_h_s_weak = mean(erp_epoch_h_s_weak([channel1,channel2,channel3,channel4],:),1);

mean_s_h_weak = mean(erp_epoch_s_h_weak([channel1,channel2,channel3,channel4],:),1);
mean_s_n_weak = mean(erp_epoch_s_n_weak([channel1,channel2,channel3,channel4],:),1);
mean_s_s_weak = mean(erp_epoch_s_s_weak([channel1,channel2,channel3,channel4],:),1);

mean_n_h_weak = mean(erp_epoch_n_h_weak([channel1,channel2,channel3,channel4],:),1);
mean_n_n_weak = mean(erp_epoch_n_n_weak([channel1,channel2,channel3,channel4],:),1);
mean_n_s_weak = mean(erp_epoch_n_s_weak([channel1,channel2,channel3,channel4],:),1);

% Unbewusst
mean_h_h_strong = mean(erp_epoch_h_h_strong([channel1,channel2,channel3,channel4],:),1);
mean_h_n_strong = mean(erp_epoch_h_n_strong([channel1,channel2,channel3,channel4],:),1);
mean_h_s_strong = mean(erp_epoch_h_s_strong([channel1,channel2,channel3,channel4],:),1);

mean_s_h_strong = mean(erp_epoch_s_h_strong([channel1,channel2,channel3,channel4],:),1);
mean_s_n_strong = mean(erp_epoch_s_n_strong([channel1,channel2,channel3,channel4],:),1);
mean_s_s_strong = mean(erp_epoch_s_s_strong([channel1,channel2,channel3,channel4],:),1);

mean_n_h_strong = mean(erp_epoch_n_h_strong([channel1,channel2,channel3,channel4],:),1);
mean_n_n_strong = mean(erp_epoch_n_n_strong([channel1,channel2,channel3,channel4],:),1);
mean_n_s_strong = mean(erp_epoch_n_s_strong([channel1,channel2,channel3,channel4],:),1);

%% Save ERP
% Bewusst
Erp_P100.h_h_weak = mean_h_h_weak;
Erp_P100.h_n_weak = mean_h_n_weak;
Erp_P100.h_s_weak = mean_h_s_weak;

Erp_P100.s_h_weak = mean_s_h_weak;
Erp_P100.s_n_weak = mean_s_n_weak;
Erp_P100.s_s_weak = mean_s_s_weak;

Erp_P100.n_h_weak = mean_n_h_weak;
Erp_P100.n_n_weak = mean_n_n_weak;
Erp_P100.n_s_weak = mean_n_s_weak;

% Unbewusst
Erp_P100.h_h_strong = mean_h_h_strong;
Erp_P100.h_n_strong = mean_h_n_strong;
Erp_P100.h_s_strong = mean_h_s_strong;

Erp_P100.s_h_strong = mean_s_h_strong;
Erp_P100.s_n_strong = mean_s_n_strong;
Erp_P100.s_s_strong = mean_s_s_strong;

Erp_P100.n_h_strong = mean_n_h_strong;
Erp_P100.n_n_strong = mean_n_n_strong;
Erp_P100.n_s_strong = mean_n_s_strong;


%% Get P100 and Save

% Bewusst
h_h_weak_P100Range = mean_h_h_weak(1,(range_min:range_max)); 
h_n_weak_P100Range = mean_h_n_weak(1,(range_min:range_max)); 
h_s_weak_P100Range = mean_h_s_weak(1,(range_min:range_max)); 

s_h_weak_P100Range = mean_s_h_weak(1,(range_min:range_max)); 
s_n_weak_P100Range = mean_s_n_weak(1,(range_min:range_max)); 
s_s_weak_P100Range = mean_s_s_weak(1,(range_min:range_max)); 

n_h_weak_P100Range = mean_n_h_weak(1,(range_min:range_max)); 
n_n_weak_P100Range = mean_n_n_weak(1,(range_min:range_max)); 
n_s_weak_P100Range = mean_n_s_weak(1,(range_min:range_max)); 

% Unbewusst
h_h_strong_P100Range = mean_h_h_strong(1,(range_min:range_max)); 
h_n_strong_P100Range = mean_h_n_strong(1,(range_min:range_max)); 
h_s_strong_P100Range = mean_h_s_strong(1,(range_min:range_max)); 

s_h_strong_P100Range = mean_s_h_strong(1,(range_min:range_max)); 
s_n_strong_P100Range = mean_s_n_strong(1,(range_min:range_max)); 
s_s_strong_P100Range = mean_s_s_strong(1,(range_min:range_max)); 

n_h_strong_P100Range = mean_n_h_strong(1,(range_min:range_max)); 
n_n_strong_P100Range = mean_n_n_strong(1,(range_min:range_max)); 
n_s_strong_P100Range = mean_n_s_strong(1,(range_min:range_max)); 

%% Save P100
% Bewusst 
Erp_P100.h_h_weak_P100 = min(h_h_weak_P100Range);
Erp_P100.h_n_weak_P100 = min(h_n_weak_P100Range);
Erp_P100.h_s_weak_P100 = min(h_s_weak_P100Range);

Erp_P100.s_h_weak_P100 = min(s_h_weak_P100Range);
Erp_P100.s_n_weak_P100 = min(s_n_weak_P100Range);
Erp_P100.s_s_weak_P100 = min(s_s_weak_P100Range);

Erp_P100.n_h_weak_P100 = min(n_h_weak_P100Range);
Erp_P100.n_n_weak_P100 = min(n_n_weak_P100Range);
Erp_P100.n_s_weak_P100 = min(n_s_weak_P100Range);

% Unbewusst 
Erp_P100.h_h_strong_P100 = min(h_h_strong_P100Range);
Erp_P100.h_n_strong_P100 = min(h_n_strong_P100Range);
Erp_P100.h_s_strong_P100 = min(h_s_strong_P100Range);

Erp_P100.s_h_strong_P100 = min(s_h_strong_P100Range);
Erp_P100.s_n_strong_P100 = min(s_n_strong_P100Range);
Erp_P100.s_s_strong_P100 = min(s_s_strong_P100Range);

Erp_P100.n_h_strong_P100 = min(n_h_strong_P100Range);
Erp_P100.n_n_strong_P100 = min(n_n_strong_P100Range);
Erp_P100.n_s_strong_P100 = min(n_s_strong_P100Range);



%% Only Primer
% Get Index of stimulus
index_h_weak = [index_h_h_weak, index_h_n_weak, index_h_s_weak];
index_n_weak = [index_n_h_weak, index_n_n_weak, index_n_s_weak];
index_s_weak = [index_s_h_weak, index_s_n_weak, index_s_s_weak];
index_h_strong = [index_h_h_strong, index_h_n_strong, index_h_s_strong];
index_n_strong = [index_n_h_strong, index_n_n_strong, index_n_s_strong];
index_s_strong = [index_s_h_strong, index_s_n_strong, index_s_s_strong];

index_h = [index_h_h_weak, index_h_n_weak, index_h_s_weak, index_h_h_strong, index_h_n_strong, index_h_s_strong];
index_n = [index_n_h_weak, index_n_n_weak, index_n_s_weak, index_n_h_strong, index_n_n_strong, index_n_s_strong];
index_s = [index_s_h_weak, index_s_n_weak, index_s_s_weak, index_s_h_strong, index_s_n_strong, index_s_s_strong];

index_weak = [index_h_h_weak, index_h_n_weak, index_h_s_weak, index_n_h_weak, index_n_n_weak, index_n_s_weak, index_s_h_weak, index_s_n_weak, index_s_s_weak];
index_strong = [index_h_h_strong, index_h_n_strong, index_h_s_strong, index_n_h_strong, index_n_n_strong, index_n_s_strong, index_s_h_strong, index_s_n_strong, index_s_s_strong];

%Get Trials with respective index
trials_h_weak = Epoch(index_h_weak);
trials_n_weak = Epoch(index_n_weak);
trials_s_weak = Epoch(index_s_weak);
trials_h_strong = Epoch(index_h_strong);
trials_n_strong = Epoch(index_n_strong);
trials_s_strong = Epoch(index_s_strong);

trials_happy = Epoch(index_h);
trials_neutral = Epoch(index_n);
trials_sad = Epoch(index_s);

trials_weak = Epoch(index_weak);
trials_strong = Epoch(index_strong);

%Create a mean ERP for all trials
erp_epoch_h_weak = mean(EEG.data(:,:,trials_h_weak),3);
erp_epoch_n_weak = mean(EEG.data(:,:,trials_n_weak),3);
erp_epoch_s_weak = mean(EEG.data(:,:,trials_s_weak),3);
erp_epoch_h_strong = mean(EEG.data(:,:,trials_h_strong),3);
erp_epoch_n_strong = mean(EEG.data(:,:,trials_n_strong),3);
erp_epoch_s_strong = mean(EEG.data(:,:,trials_s_strong),3);

erp_epoch_happy = mean(EEG.data(:,:,trials_happy),3);
erp_epoch_neutral = mean(EEG.data(:,:,trials_neutral),3);
erp_epoch_sad = mean(EEG.data(:,:,trials_sad),3);

erp_epoch_weak = mean(EEG.data(:,:,trials_weak),3);
erp_epoch_strong =mean(EEG.data(:,:,trials_strong),3);


%Create a mean for the Channels we are interested in
mean_h_weak = mean(erp_epoch_h_weak([channel1,channel2,channel3,channel4],:),1);
mean_n_weak = mean(erp_epoch_n_weak([channel1,channel2,channel3,channel4],:),1);
mean_s_weak = mean(erp_epoch_s_weak([channel1,channel2,channel3,channel4],:),1);
mean_h_strong = mean(erp_epoch_h_strong([channel1,channel2,channel3,channel4],:),1);
mean_n_strong = mean(erp_epoch_n_strong([channel1,channel2,channel3,channel4],:),1);
mean_s_strong = mean(erp_epoch_s_strong([channel1,channel2,channel3,channel4],:),1);

mean_h = mean(erp_epoch_happy([channel1,channel2,channel3,channel4],:),1);
mean_n = mean(erp_epoch_neutral([channel1,channel2,channel3,channel4],:),1);
mean_s = mean(erp_epoch_sad([channel1,channel2,channel3,channel4],:),1);

mean_weak = mean(erp_epoch_weak([channel1,channel2,channel3,channel4],:),1);
mean_strong = mean(erp_epoch_strong([channel1,channel2,channel3,channel4],:),1);

% Save Erp
Erp_P100.h_weak = mean_h_weak;
Erp_P100.n_weak = mean_n_weak;
Erp_P100.s_weak = mean_s_weak;
Erp_P100.h_strong = mean_h_strong;
Erp_P100.n_strong = mean_n_strong;
Erp_P100.s_strong = mean_s_strong;

Erp_P100.happy = mean_h;
Erp_P100.neutral = mean_n;
Erp_P100.sad = mean_s;

Erp_P100.weak = mean_weak;
Erp_P100.strong = mean_strong;
        
    % Effects
Erp_P100.cons_effect = mean_weak-mean_strong;
Erp_P100.cons_effect_happy = mean_h_weak-mean_h_strong;
Erp_P100.cons_effect_neutral = mean_n_weak-mean_n_strong;
Erp_P100.cons_effect_sad = mean_s_weak-mean_s_strong;

Erp_P100.emotion_effect_happy_conscious =  mean_n_weak-mean_h_weak;
Erp_P100.emotion_effect_sad_conscious = mean_n_weak-mean_s_weak;
Erp_P100.emotion_effect_happy_unconscious = mean_n_strong-mean_h_strong;
Erp_P100.emotion_effect_sad_unconscious = mean_n_strong-mean_s_strong;




% Get range of P100
h_weak_P100Range = mean_h_weak(1,(range_min:range_max));
n_weak_P100Range = mean_n_weak(1,(range_min:range_max)); 
s_weak_P100Range = mean_s_weak(1,(range_min:range_max)); 
h_strong_P100Range = mean_h_strong(1,(range_min:range_max));
n_strong_P100Range = mean_n_strong(1,(range_min:range_max)); 
s_strong_P100Range = mean_s_strong(1,(range_min:range_max)); 


happy_P100Range = mean_h(1,(range_min:range_max));
neutral_P100Range = mean_n(1,(range_min:range_max)); 
sad_P100Range = mean_s(1,(range_min:range_max)); 

weak_P100Range = mean_weak(1,(range_min:range_max)); 
strong_P100Range = mean_strong(1,(range_min:range_max));

cons_effect_P100Range = Erp_P100.cons_effect(1,(range_min:range_max)); 
cons_effect_happy_P100Range = Erp_P100.cons_effect_happy (1,(range_min:range_max)); 
cons_effect_neutral_P100Range = Erp_P100.cons_effect_neutral(1,(range_min:range_max)); 
cons_effect_sad_P100Range = Erp_P100.cons_effect_sad(1,(range_min:range_max)); 

emotion_effect_happy_conscious_P100Range =  Erp_P100.emotion_effect_happy_conscious(1,(range_min:range_max)); 
emotion_effect_sad_conscious_P100Range = Erp_P100.emotion_effect_sad_conscious(1,(range_min:range_max)); 
emotion_effect_happy_unconscious_P100Range = Erp_P100.emotion_effect_happy_unconscious(1,(range_min:range_max)); 
emotion_effect_sad_unconscious_P100Range = Erp_P100.emotion_effect_sad_unconscious(1,(range_min:range_max)); 

% Save P100
Erp_P100.h_weak_P100 = min(h_weak_P100Range);
Erp_P100.n_weak_P100 = min(n_weak_P100Range);
Erp_P100.s_weak_P100 = min(s_weak_P100Range);
Erp_P100.h_strong_P100 = min(h_strong_P100Range);
Erp_P100.n_strong_P100 = min(n_strong_P100Range);
Erp_P100.s_strong_P100 = min(s_strong_P100Range);

Erp_P100.happy_P100= min(happy_P100Range);
Erp_P100.neutral_P100 = min(neutral_P100Range);
Erp_P100.sad_P100 = min(sad_P100Range);

Erp_P100.weak_P100 = min(weak_P100Range);
Erp_P100.strong_P100 = min(strong_P100Range);

    % Effects
Erp_P100.cons_effect_P100 = min(cons_effect_P100Range);
Erp_P100.cons_effect_happy_P100 = min(cons_effect_happy_P100Range);
Erp_P100.cons_effect_sad_P100 = min(cons_effect_neutral_P100Range);
Erp_P100.cons_effect_neutral_P100 = min(cons_effect_sad_P100Range);

Erp_P100.emotion_effect_happy_conscious_P100 = min(emotion_effect_happy_conscious_P100Range);
Erp_P100.emotion_effect_sad_conscious_P100 = min(emotion_effect_sad_conscious_P100Range);
Erp_P100.emotion_effect_happy_unconscious_P100 = min(emotion_effect_happy_unconscious_P100Range);
Erp_P100.emotion_effect_sad_unconscious_P100 = min(emotion_effect_sad_unconscious_P100Range);

% Zum Plotten
Erp_P100.times = EEG.times;



