# Unconscious Bias in Major Depressive Disorder: EEG-fMRI Study

[![DOI](https://img.shields.io/badge/DOI-10.1016/j.bpsc.2024.07.005-blue)](https://doi.org/10.1016/j.bpsc.2024.07.005)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC--BY--4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![Status: Open Access](https://img.shields.io/badge/Status-Open%20Access-brightgreen.svg)](https://doi.org/10.17605/OSF.IO/37XD2)

---

#### AUTHORS 
Julia Schräder 1,2, Lennard Herzberg 1, Han-Gue Jo 3, Lucia Hernandez-Pena 1,2, Julia Koch 1,2, Ute Habel 1,2, Lisa Wagels 1,2 

* <sub><sup>1 Department of Psychiatry, Psychotherapy and Psychosomatics, Medical Faculty, Uniklinik RWTH Aa-chen University, Pauwelstraße 30, 52074 Aachen, Germany</sup></sub>
*  <sub><sup>2 Institute of Neuroscience and Medicine: JARA-Institute Brain Structure Function Relationship (INM 10), Research Center Jülich, Jülich, Germany</sup></sub>
*  <sub><sup>3 School of Computer SoftwareInformation and Communication Engineering, Kunsan National University, 588 Daehak-ro Gunsan, South Korea</sup></sub>

---

## General Information

This study investigates the perception of emotional stimuli presented at an unconscious level in a sample of 66 healthy individuals and 60 patients with major depressive disorder (MDD).  

To examine the neurophysiological pathways of biased unconscious emotion processing in depression, we applied **simultaneous EEG-fMRI** recordings, focusing on the event-related potential for facial expressions (*N170*) to predict the presence of unconscious cognitive biases in depressed individuals.  

A **modified backward mask priming task** was used, presenting facial expressions (happy, sad, neutral) under both unconscious (16.7 ms) and conscious (150 ms) prime durations, with congruent or incongruent target emotions.  

---
## ✨ Abstract


**Background:** Major depressive disorder (MDD) is characterized by strong emotional dysregulation. Mechanisms driving the negative affect in depression may be fast processes existing on an unconscious level.
**Methods:** A priming task was conducted using simultaneous electroencephalography–functional magnetic resonance imaging measurement involving presentation of facial expressions (happy, sad, and neutral) to examine the neurophysiological pathway of biased unconscious emotion processing in MDD. Priming prior to a target emotion created unconscious (16.7-ms primer) and conscious (150-ms primer) trials. A large sample (N = 126) was recruited, containing healthy control participants (n = 66; 37 women) and participants with MDD (n = 60; 31 women).
**Results:** The healthy control group showed a shorter reaction time in happy but not in sad or neutral trials compared with the MDD group. N170 amplitudes were lower in trials with unconscious than conscious primer presentation. N170 amplitudes correlated with cortical (right fusiform gyrus, right middle temporal gyrus, right inferior temporal gyrus, left supplementary motor area, right middle frontal gyrus) and subcortical brain regions (right amygdala). The strength of N170 and brain activity correlation increased when the stimulus was consciously presented. Presented emotions did not affect the correlation of N170 values and brain activity.
**Conclusions:** Our findings show that MDD may exhibit biased emotion regulation abilities at a behavioral and neurophysiological level. Face-sensitive event-related potentials demonstrate a correlation with heightened brain activity in regions associated with both face recognition (fusiform gyrus) and emotion processing (amygdala). These findings are evident in both MDD and healthy control groups, with lower effect sizes in the MDD group indicating reduced emotion recognition and processing abilities.

---

## Repository Contents

- `data/` – datasets for GLMM calculation, questionnaire results, and mean reaction time data  
- `scripts/` – analysis scripts and datasets for GLMM, EEG, fMRI, and EEG-informed fMRI analysis  

---

## Experimental Design

The backward-masked priming paradigm included 36 gender-balanced images (12 happy, 12 sad, 12 neutral) from the **FACES database**. Each image was shown 10 times on a 120 Hz LCD monitor using **PsychoPy 3**.  

- **Fixation cross**: 300 ms (36 frames)  
- **Prime**:  
  - Strongly masked: 16.7 ms prime + 66.7 ms mask  
  - Weakly masked: 150 ms prime + 66.7 ms mask  
- **Target**: emotional facial expression (happy, sad, neutral)  
- **Response phase**: 1.5 s; participants responded via button box:  
  - Index finger = sad  
  - Middle finger = neutral  
  - Ring finger = happy  

This paradigm ensured distinct **unconscious** and **conscious** emotional conflict vs. non-conflict trials to test unconscious emotional processing effects.  




---

## Pre-Registration

This project is openly available via OSF:  
- OSF Project Page: [https://osf.io/37xd2](https://osf.io/37xd2)  
- DOI: [10.17605/OSF.IO/37XD2](https://doi.org/10.17605/OSF.IO/37XD2)  

---

## Technologies Used

- PsychoPy3 v2020.2.4  
- RStudio v2023.06.1  
- Python 3  
- GraphPad Prism v10.0.2  
- MATLAB R2020b  
- EEGLAB v2022.1  

---

## Acknowledgements

This study was supported by the **Brain Imaging Facility of the Interdisciplinary Center for Clinical Research**, Faculty of Medicine, RWTH Aachen University, and the **International Research Training Group (IRTG 2150)**.  

Funding:  
- Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) — 667892 / JO 1453/2-1  
- International Research Training Group (IRTG 2150) — 269953372 / GRK2150  
- FZJ-NST Bilateral Cooperation Program funded by Forschungszentrum Jülich and the National Research Council of Science & Technology (Global-22-001)  

---

## 📜 License

This work is distributed under the terms of the  
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).  
