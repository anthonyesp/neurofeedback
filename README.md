# Visual and haptic feedback in detecting motor imagery within a wearable brain-computer interface

## Description
This repository contains data and code about the detection of motor imagery within a wearable brain-computer interface when visual and haptic feedbacks were provided to the user.

Eight right-handed volunteers (three males, mean age 28 years, and five females, mean age 25 years) participated in the experiments. These were conducted at the Augmented Reality for Health Monitoring Laboratory (ARHeMLab, University of Naples Federico II) in Italy.

The experiments were carried out in two or three sessions on different days, each lasting about two hours. Participants were asked to imagine the movement of the left or right hand.

In each session, participants first performed motor imagery without any feedback (B0 block), which was required for identifying the online classification model. Then, they received online feedback (B1_h, B1_v, and B1_vh for haptic, visual and visual-haptic feedback respectively). The presentation of the three feedback modalities were randomized for each subject to avoid biases associated with the sequence of presentation.

The blocks consisted of three runs with 30 trials each (15 per class) with about two-minute breaks in between. The order of the cue-based motor imagery was again randomized to avoid any bias. The timing of a single trial was recalled from the standard paradigms of BCI competitions. In particular, it consists of an initial relax, a cue at t = 2 s indicating the task to carry on, motor imagery starting at t = 3 s, and motor imagery ending at t = 6 s. Final relaxation was then presented, and its duration was randomized between 1 s to 2 s.

Block B0 could eventually contain a fourth run if the classification accuracy associated with the first three runs was not satisfying.

The scripts provided here were implemented in Matlab R2021a, but they could be used with previous versions of Matlab as well. Data processing relies on the "filter bank common spatial pattern" algorithm, and inherent scripts must be firstly downloaded here: https://github.com/anthonyesp/channel_selection.git (only subfolder “FBCSP algorithm”).

## Data file description
Data are stored with the .mat extension. For each session, block, and subject, a folder contains three or four .mat files related with the recorded runs. Once the .mat file is loaded, the workspace will contain:

•	for block B0: a cell “data”, a “nf” integer value representing a counter for the experimental runs, and a structure “notes” with meta-information.

•	for neurofeedback blocks: a cell “data”, a structure “notes” with meta-information, a structure “tr_param” with the parameters of the model used in online classification (obtained with offline training on B0 data), and a char “type” with the information related to the feedback type.

Each “data” contains: 

–	X : the EEG data stream given as a M x N matrix, where N is the number of channels and M is the number of samples; 

–	trial : a vector of length 30 containing the starting sample per each trial; 

–	y: a vector of length 30 containing the labels per each trial; 

–	fs: the sampling frequency; 

–	classes: a 1 × 2 cell containing the motor imagery task associated to the two possible labels (1 and 2); 

–	artifacts: a logical vector of length 30. It flags artifact marked by visual inspection. Specifically, 0 corresponds to a clean trial and 1 corresponds to a trial containing an artifact.

–	gender of the participant.

–	age of the participant. 

## Furter details
1. Arpaia, P., Coyle, D., Donnarumma, F., Esposito, A., Natalizio, A., & Parvis, M. (2023). Visual and haptic feedback in detecting motor imagery within a wearable brain–computer interface. Measurement, 206, 112304.

2. "Motor imagery brain-computer interface and extended reality" project on ResearchGate
