# Unconscious Bias in Major Depressive Disorder: EEG-fMRI Study

## Table of contents
* [General info](#general-info)
	* [Authors](#authors) 
	* [Publication](#OSF) 
	* [License](#license)
* [Experimental design](#experimental-design)
* [Folder descripion](#folder-description)
* [Technologies](#technologies)


## General info
This study is about the perception of emotional stimuli presented at an unconscious level in a sample of 66 healthy individuals and patients 60 with major depressive disorder. To examine the neurophysiological pathway of biased unconscious emotion processing in depression, simultaneous EEG-fMRI measures are applied focusing on the event-related potential for facial expressions (N170) for predicting the presence of unconscious cognitive biases in depressed individuals.
A modified backward mask priming task was conducted utilizing simultaneous EEG-fRMI measurement involving presentation of facial expressions (happy, sad, neutral). Priming prior to a target emotion created distinct emotional trials: unconscious (16.7 ms primer duration) and conscious (150 ms primer duration) with either congruency or incongruency of emotions between primer and target. 


#### AUTHORS
Julia Schräder 1,2, Lennard Herzberg 1, Han-Gue Jo 3, Lucia Hernandez-Pena 1,2, Julia Koch 1,2,, Ute Habel 1,2, Lisa Wagels 1,2

* <sub><sup>1 Department of Psychiatry, Psychotherapy and Psychosomatics, Medical Faculty, Uniklinik RWTH Aa-chen University, Pauwelstraße 30, 52074 Aachen, Germany</sup></sub>
* <sub><sup>2 Institute of Neuroscience and Medicine: JARA-Institute Brain Structure Function Relationship (INM 10), Research Center Jülich, Jülich, Germany</sup></sub>
* <sub><sup>3 School of Computer Information and Communication Engineering, Kunsan National University, Gun-san, Korea</sup></sub>


#### OSF 

[osf.io/mscz5](https://osf.io/37xd2)

#### License

CC-By Attribution 4.0 International 

## Experimental design

The backward mask priming paradigm provided masked stimulus presentation of emotional and neutral stimuli. Images of 12 happy, 12 sad and 12 neutral facial expressions were presented. A total of 36 gender balanced images taken from the FACES database were used. Each image was presented 10 times against a grey background at the center of an LCD monitor (screen refresh rate = 120 Hz) using the PsychoPy 3 software.

Every trial started with a fixation cross presented for 300 ms (36 frames). Afterwards, a prime stimulus appeared either strongly masked (with a presentation time of 16.7 ms followed by a mask for 66.7 ms) or weakly masked (with a primer for 150 ms followed by a mask for 66.7 ms). Participants were asked to rate the emotional expression of the target image via a button box using the index finger to respond “sad”, the middle finger to respond “neutral”, or the ring finger to respond “happy”. Participants had a response phase of 1.5 s. 

Using either images of the same (no-conflict) or different (conflict) emotional expression (happy, sad, neutral) as primer / target, this modified backward-masked priming task ensured two distinct unconscious (16.7 ms) and conscious (150 ms) emotional conflict trials as well as non-conflict trials to test the effect of unconscious emotional information on cognitive processing.

![Task](https://github.com/JuliaSchraeder/UnconsciousBias/assets/54576554/270e406b-83e8-4430-ac6d-9df84326ba85)

## Folder description

* `data` includes data used for GLMM calculation, questionnaire results and mean reaction time data. 
* `scripts` includes script and dataset for GLMM, EEG, fMRI and eeg-informed fMRI anaylsis 

## Technologies
Project is created with:
* PsychoPy3: Version v2020.2.4
* RStudio: version 2023.06.1
* Python 3
* GraphPad Prism 10.0.2.
* Matlab version 2020b


