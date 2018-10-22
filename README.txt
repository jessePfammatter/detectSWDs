README for ‘detectSWD’ (SWD = spike and wave discharge) set of scripts. (updated 08/22/18)

This repository is associated with the manuscript, 'An automated, machine learning-based detection algorithm for spike-wave discharges (SWDs) in a mouse model of absence epilepsy' which can be found at https://doi.org/10.1101/309146. For questions please contact Jesse Pfammatter at jesse.pfammatter@gmail.com or pfammatter@wisc.edu

% ----- Demo Materials ----- %

% ----- Training a Model ----- %


1. A demo

1. Obtain EEG from mouse
2. Get some files tagged (4 second epochs) for SWDs in Svarog or Serenia Sleep Pro to help train the model.
3. open the developSWDClassifier.m script and modify it to include the files that you have tagged for SWDs. Use this script to build a Support Vector Machine (SVM) model that can identify SWDs in a new set of files. There are some test features built into the script that will allow you to see the number of true positives, false positives, etc. When developing a model, shoot for the smallest number of false negatives as the automated portion of the script will go through and finish the job and properly reduce the false positives (also, there might be more SWD activity than you think). An initial model that missed a lot of true positives will not be able to find them during the exact detection phase (as it only looks at epochs flagged by the SVM).

We have used this script to identify SWDs in mice with the R43Q mutation which causes absence epilepsy. These models were developed from a series of 24 hour EEG files that were
collected by Aaron Nelson and Kile Mangon ~ 2013.

There is also a set of files in characterizeSWDs that were used to develop the SWD template and will be modified to further classify types of SWDs that detectSWDs_automated outputs.

This algorithm is advantageous because it allows a framework for scoring SWDs that removes
bias based on sleep staging, reduces errors due to variation in scorer attentiveness and
skill level and, reduces the time needed to score a file. 

This analysis packages was developed by Jesse Pfammatter and Mathew Jones (2016).

% ----- What's in this folder?

-------
*detectSWDs_automated.m*  
** START HERE IF YOU ALREADY HAVE A TRAINED CLASSIFIER
This is the main script used to analyze new EEG records for SWDs and is all you need to develop a list of exact SWD locations from an EEG file if you already have a trained classifier.
-------

-------
*developSWDClassifier.m*  
** START HERE IF YOU DON'T HAVE A CLASSIFIER TRAINED
This script is a wrapper for all the code needed to develop a new SVM. It's where to start if you don't already have a trained classifier. It’s not super user friendly at the current moment, so this needs to be improved if we end up incorporating it into the GUI. At the current moment it does not need to be placed into a GUI.
-------

-------
*generateSWDPredictors.m*
This script takes inputs for EDF and seizure scorings. See the help file for the script for more information.
-------

-------
* trainPolynomialSVMClassifier.m
This is the code that trains the actual SVM and where tweaks to the model parameters should be made if needed.
-------

-------
SWD detection GUI folder
Matt was working on this script, and might continue at some point if we get that grant.
-------

-------
Auxiliary scripts developed initially for this project, but that eventually were moved to
a utilities folder or something (or will be if they have not already). This list is not complete yet.

* normalizeEEG.m
* specifyExactSWDLocations.m
-------
