% Test detectSWDs

%% variable argument defaults

% [basePath, ~, ~] = getUserPath();
[basePath, ~, cookieMonster] = getUserPath();

% edffilespec = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/EDFs/Ronde_4.edf');
edffilespec = '/Users/jonesmat/Desktop/TrashMe/TempData/Ronde_4.edf'


%% load EDF file
EDF = read_EDF_mj(edffilespec);


%% calculate predictive matrix
disp('Generating the Matrix of Predictors...');
predictorMatrix = generateSWDPredictors_mj( EDF );


%% load SVM classifier
% savepath = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/SWDClassificationData/');
savepath = '/Users/Shared/JonesLab_Git/joneslab_code/Projects/sleep_and_seizures/EEG_Analysis/seizureIdentification/SVMClassifier_SWDs/detectSWDs GUI/'
load(strcat(savepath, 'polynomialSWDSVMClassifier.mat'));

%% classify the SWDs
[yfit, ~] = trainedClassifier.predictFcn(predictorMatrix); 
tmpfig = figure;
plot(yfit)




%% Detect SWDs
% detectSWDs(edffilespec);
% detectSWDs_mj( edffilespec);
