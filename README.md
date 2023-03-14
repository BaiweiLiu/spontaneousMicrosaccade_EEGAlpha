# spontaneousMicrosaccade_EEGAlpha
Code associated with the article "Microsaccades transiently lateralise EEG alpha activity" by Baiwei Liu, Anna C. Nobre, Freek van Ede. 

date: 14 Mar 2023. 

contact: baiwei liu (b.liu@vu.nl)

These codes are written in Matlab R2021.

To run these codes, you can change the path where you put the code and the data at the beginning of these codes.

You can run these codes in following order:
1. Install necessary toolbox: fieldtrip and gramm and then install our custom toolbox: futureToolBox
2. run eeg_preprocess.m and eye_preprocess.m to preprocess eeg and eye data 
3. run cloudPipPO78O12.m to get all the results used in the paper. 
4. run the codes in code4figures to get the figure in the paper.

The behavioural and eye-tracking data can be freely downloaded at: https://doi.org/10.5061/dryad.m99r286 (Experiment 1). 
The corresponding EEG data can be freely downloaded at: https://doi.org/10.5061/dryad.sk8rb66.
