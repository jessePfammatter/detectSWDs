
% ----- import the namefile and load the aggregated data

namefilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleep_and_seizures_RQ_namefile.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'
savePath = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/seizureScoring/SWDtest3/');
load(strcat(savePath, 'SWDDATA_07_27_17_normalized.mat')); 

% ----- generate event specific predictors and assemble PCA matrix including all variables, even ones that don't calculate properly if no peaks are found -- those we will have to drop from the analysis
fullMat = generatePredictorsForEventClassifier(data, 256);

% ----- pull responses from data
% Rama's full set of scoring
load(strcat(savePath, 'data/SWDDATA_07_27_17_RM.mat'));
%load(strcat(savePath, 'data/SWDGUIDATA_4-24-17_RM.mat')); % data from SWDtest1
dataS1 = D.data;
for i = 1:size(data, 2)
    labelsS1(i) = dataS1(i).responses;
end

% Matt's full set of scoring
load(strcat(savePath, 'data/SWDDATA_07_27_17_MJ.mat'));
%load(strcat(savePath, 'data/SWDGUIDATA_4-24-17_MJ.mat')); % data from SWDtest1
dataS2 = D.data;
for i = 1:size(data, 2)
    labelsS2(i) = dataS2(i).responses;
end

% Kile's responses
load(strcat(savePath, 'data/SWDDATA_07_27_17_KPM.mat')); 
dataS3 = D.data;
for i = 1:size(data, 2)
    labelsS3(i) = dataS3(i).responses;
end

% Jesse's responses
load(strcat(savePath, 'data/SWDDATA_07_27_17_JP.mat')); 
%load(strcat(savePath, 'data/SWDGUIDATA_4-24-17_JP.mat')); % data from SWDtest1
dataS4 = D.data;
for i = 1:size(data, 2)
    labelsS4(i) = dataS4(i).responses;
end

% make a vector of the joint responses
labelsHC = [];
for i = 1:size(data, 2)
    if labelsS1(i) == 2 & labelsS2(i) == 2 & labelsS3(i) == 2 & labelsS4(i) == 2
        labelsHC(i) = 2; %{'SWD (Agreed)'};  
    elseif labelsS1(i) == 1 & labelsS2(i) == 1 & labelsS3(i) == 1 & labelsS4(i) == 1
        labelsHC(i) = 1; %{'nonSWD (Agreed)'};
    else
        labelsHC(i) = 0; %{'disagreed'};
    end
end

% make color labels for each of labelsS1-3 and for labelsHCs
for i = 1:length(labelsS1)
    if labelsS1(i) == 2
        labelsS1Colors(i,:) = [0, 0, 255]; % cyan rgb
    elseif labelsS1(i) == 1
        labelsS1Colors(i,:) = [255, 0, 0]; % red rgb
    else
        labelsS1Colors(i,:) = [255, 255, 0];
    end
end

for i = 1:length(labelsS2)
    if labelsS2(i) == 2
        labelsS2Colors(i,:) = [0, 0, 255]; % cyan rgb
    elseif labelsS2(i) == 1
        labelsS2Colors(i,:) = [255, 0, 0]; % red rgb
    else
        labelsS2Colors(i,:) = [255, 255, 0];
    end
end

for i = 1:length(labelsS3)
    if labelsS3(i) == 2
        labelsS3Colors(i,:) = [0, 0, 255]; % cyan rgb
    elseif labelsS3(i) == 1
        labelsS3Colors(i,:) = [255, 0, 0]; % red rgb
    else
        labelsS3Colors(i,:) = [255, 255, 0];
    end
end

for i = 1:length(labelsS4)
    if labelsS4(i) == 2
        labelsS4Colors(i,:) = [0, 0, 255]; % cyan rgb
    elseif labelsS4(i) == 1
        labelsS4Colors(i,:) = [255, 0, 0]; % red rgb
    else
        labelsS4Colors(i,:) = [255, 255, 0];
    end
end

for i = 1:length(labelsHC)
    if labelsHC(i) == 2
        labelsHCColors(i,:) = [0, 0, 255]; % cyan rgb
    elseif labelsHC(i) == 1
        labelsHCColors(i,:) = [255, 0, 0]; % red rgb
    else
        labelsHCColors(i,:) = [255, 255, 0];
    end
end

% create the weighted labels data
clear weightedResponseLabels weightedResponseWeights additionalIndex
for i = 1:length(labelsHC)
     
    if labelsS1(i) == 0
        templabelsS1(i) = 1.5;
    else
        templabelsS1(i) = labelsS1(i);
    end
    if labelsS2(i) == 0
        templabelsS2(i) = 1.5;
    else
        templabelsS2(i) = labelsS2(i);
    end    
    if labelsS3(i) == 0
        templabelsS3(i) = 1.5;
    else
        templabelsS3(i) = labelsS3(i);
    end
    if labelsS4(i) == 0
        templabelsS4(i) = 1.5;
    else
        templabelsS4(i) = labelsS4(i);
    end     
    if mean([templabelsS1(i), templabelsS2(i), templabelsS3(i), templabelsS4(i)]) > 1.5
        weightedResponseLabels(i) = 2;
        if templabelsS1(i) == 2
            s1w = 1;
        else
            s1w = 0;
        end
        if templabelsS2(i) == 2
            s2w = 1;
        else
            s2w = 0;
        end
        if templabelsS3(i) == 2
            s3w = 1;
        else
            s3w = 0;
        end
        if templabelsS4(i) == 2
            s4w = 1;
        else
            s4w = 0;
        end
        weightedResponseWeights(i) = (s1w + s2w + s3w + s4w) / 4;
    else
        weightedResponseLabels(i) = 1;
        if templabelsS1(i) == 1
            s1w = 1;
        else
            s1w = 0;
        end
        if templabelsS2(i) == 1
            s2w = 1;
        else
            s2w = 0;
        end
        if templabelsS3(i) == 1
            s3w = 1;
        else
            s3w = 0;
        end
        if templabelsS4(i) == 1
            s4w = 1;
        else
            s4w = 0;
        end
        weightedResponseWeights(i) = (s1w + s2w + s3w + s4w) / 4;
    end
end
 
% make color labels for weighted labels
for i = 1:length(weightedResponseLabels)
    if weightedResponseLabels(i) == 2
        weightedResponseLabelsColors(i,:) = [0, 0, 255]; % cyan rgb
    elseif weightedResponseLabels(i) == 1
        weightedResponseLabelsColors(i,:) = [255, 0, 0]; % red rgb
    end
end
 
for i = 1:length(weightedResponseLabels)
    if weightedResponseLabels(i) == 2
        weightedResponseLabelsAltColors(i,:) = [0, 0, 255]; % blue rgb
    elseif weightedResponseLabels(i) == 1
        weightedResponseLabelsAltColors(i,:) = [255, 255, 0]; % yellow rgb
    end
end
 
% clean up a few objects
clear dataS1 dataS2 dataS3 dataS4 D.data

% find variables that have NA or INF
rmNAs = find(sum(~isnan(table2array(fullMat))) == size(fullMat, 1));
rmInf = find(sum(~isinf(table2array(fullMat))) == size(fullMat, 1));
rmNAInf = intersect(rmNAs, rmInf);
fullMat = fullMat(:,rmNAInf);

% minimal training matrix with only the human concensus events
keepInd = find(labelsHC ~= 0);
HCTrainingMat = fullMat(keepInd, :);
HCTrainingLabels = labelsHC(keepInd);
HCTrainingMat.HCTrainingLabels = HCTrainingLabels';

variablesvariables = 1:size(table2array(fullMat), 2); % [1:3, 7:9];

% create a giant data matrix will all events and all responses...
humanMat = [fullMat; fullMat; fullMat; fullMat];
humanMat.labels = [labelsS1, labelsS2, labelsS3, labelsS4]';

%%
% train model
[trainedClassifier, validationAccuracy_mean, validationAccuracy_std] = trainSWDClassifier_gaussian(humanMat, 1.5, 10, 10);
%[trainedClassifier, validationAccuracy_mean, validationAccuracy_std] = trainSWDClassifier_linear(humanMat, 1.5, 1);
display(strcat({'Valdiation Accuracy with Training Data: '}, num2str(validationAccuracy_mean), {'+/-'}, num2str(validationAccuracy_std))); 

% now predict on the full data set
[labelsML, classificationScore] = trainedClassifier.predictFcn(fullMat);
classificationScore = classificationScore(:,1); % the second colum is just the inverse of the first column when you use a gaussian kernel CHANGE IF CHANGING KERNEL
labelsML = labelsML';

% make color labels for labelsML and establish cal
for i = 1:length(labelsML)
    if labelsML(i) == 2
     	labelsMLColors(i,:) = [0, 0, 255]; % cyan rgb
    else
      	labelsMLColors(i,:) = [255, 0, 0]; % red rgb
  	end
end

% add this back to the data object for the inter-rater stuff later
for i = 1:size(data, 2)
    data(i).labelsML = labelsML(i);
    data(i).classificationScore = classificationScore(i);
end

% assess quality of fit... and in relation to other scoreres
labelslist = {'S1', 'S2', 'S3', 'S4', 'ML', 'HC', 'CC'};

% add kmeans labels
labelsKM = kmeans(table2array(fullMat), 2);
labelsKM = labelsKM';

% build the table variables
fullLength = length(labelsS1);
S1_S2 = sum(labelsS1 == labelsS2) / fullLength;
S1_S3 = sum(labelsS1 == labelsS3) / fullLength;
S1_S4 = sum(labelsS1 == labelsS4) / fullLength;
S1_ML = sum(labelsS1 == labelsML) / fullLength;
S1_HC = sum(labelsS1(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
S1_CC = 1;
S2_S3 = sum(labelsS2 == labelsS3) / fullLength;
S2_S4 = sum(labelsS2 == labelsS4) / fullLength;
S2_ML = sum(labelsS2 == labelsML) / fullLength;
S2_HC = sum(labelsS2(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
S2_CC = 1;
S3_S4 = sum(labelsS3 == labelsS4) / fullLength;
S3_ML = sum(labelsS3 == labelsML) / fullLength;
S3_HC = sum(labelsS3(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
S3_CC = 1;
S4_ML = sum(labelsS4 == labelsML) / fullLength;
S4_HC = sum(labelsS4(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
S4_CC = 1;
ML_HC = sum(labelsML(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
ML_CC = 1;
HC_CC = 1;

S1_KM = sum(labelsS1 == labelsKM) / fullLength;
S2_KM = sum(labelsS2 == labelsKM) / fullLength;
S3_KM = sum(labelsS3 == labelsKM) / fullLength;
S4_KM = sum(labelsS4 == labelsKM) / fullLength;
KM_ML = sum(labelsML == labelsKM) / fullLength;
KM_HC = sum(labelsKM(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
KM_CC = 1;

matrixmatrix = [1,  S1_S2,  S1_S3,  S1_S4,  S1_KM,  S1_ML,  S1_HC,  S1_CC; ...
                0,  1,      S2_S3,  S2_S4,  S2_KM,  S2_ML,  S2_HC,  S2_CC; ...
                0,  0,      1,      S3_S4,  S3_KM,  S3_ML,  S3_HC,  S3_CC; ...
                0,  0,      0,      1,      S4_KM,  S4_ML,  S4_HC,  S4_CC; ...
                0,  0,      0,      0,      1,      KM_ML,  KM_HC,  KM_CC; ...
                0,  0,      0,      0,      0,      1,      ML_HC,  ML_CC; ...
                0,  0,      0,      0,      0,      0,      1,      HC_CC; ...
                0,  0,      0,      0,      0,      0,      0,      1];

sTable = array2table(matrixmatrix,'RowNames',labelslist,'VariableNames',labelslist) 

% CALCULATE SOME RATIOS AND SHIT HERE
averageHuman = mean([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
stdHuman = std([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
display(strcat({'Average Agreement Between Humans: '}, num2str(averageHuman), {'+/-'}, num2str(stdHuman)))
averageML_Human = mean([S1_ML, S2_ML, S3_ML, S4_ML]);
stdML_Human = std([S1_ML, S2_ML, S3_ML, S4_ML]);
display(strcat({'Average Agreement Between ML and Humans: '}, num2str(averageML_Human), {'+/-'}, num2str(stdML_Human)))

% calculate the sensitivity (true positive rate) -- should be based on the HC data
matrixmatrix = confusionmat(labelsHC, labelsML);
nEventsInHC = sum(sum(matrixmatrix(2:end, 2:end)));
specificity = matrixmatrix(2, 2) / sum(labelsHC == 1);
sensitivity = matrixmatrix(3, 3) / sum(labelsHC == 2);

%% lasso stuff after we generated the labels to see which predictors were important.

predictorNames = fullMat.Properties.VariableNames(variablesvariables);
[B,FitInfo] = lasso(table2array(fullMat), classificationScore, 'CV', 10, 'PredictorNames', predictorNames);
lassoPlot(B, FitInfo, 'PlotType', 'CV')
figure
lassoPlot(B, FitInfo)

B(:,FitInfo.Index1SE)

figure
bar(ans)
xlabel('')
%%              
% INTRA-RATER RELIABILITY FROM REPEATED EVENTS SETUP

% find all the ones that were repeated.. Should be 50 events repeated 10 times each
for i = 1:size(data, 2)
    test = data(i).signalClips(1:2048);
    for j = 1:size(data, 2)
        compare = data(j).signalClips(1:2048);
        if sum(test == compare) == 2048
            same(i, j) = 1;
        else
            same(i, j) = 0;
        end
    end
end

nReps = 10;
sameInd = find(sum(same) == nReps); % 10 total replicates

% put the responses in the data set and pull out all the events that are repeated.
for i = 1:size(data, 2)
    data(i).labelsS1 = labelsS1(i);
    data(i).labelsS2 = labelsS2(i);
    data(i).labelsS3 = labelsS3(i); 
    data(i).labelsS4 = labelsS4(i);
end
dataSame = data(sameInd);

% find the index for repeats
for i = 1:size(dataSame, 2)
    test = dataSame(i).originalSeizeNumber;
    for j = 1:size(dataSame, 2)
        compare = dataSame(j).originalSeizeNumber;
        if test == compare
            same2(i, j) = 1;
        else
            same2(i, j) = 0;
        end
    end
end

% create the events matrix
for i = 1:size(dataSame, 2)
    eventsMat(i,:) = find(same2(i,:));
end
%eventsMat

events = unique(eventsMat, 'rows');

fs = 256;
interRaterMat = generatePredictorsForEventClassifier(dataSame, fs);

[labelsMLIRM classificationScoreIRM] = trainedClassifier.predictFcn(interRaterMat);
classificationScoreIRM = classificationScoreIRM(:,1);


lengthClip = 2048;
buffer = 10;
 
for i = 1:size(dataSame, 2)
    samelabelsS1(i) = dataSame(i).labelsS1;
    samelabelsS2(i) = dataSame(i).labelsS2;
    samelabelsS3(i) = dataSame(i).labelsS3;
    samelabelsS4(i) = dataSame(i).labelsS4;
    dataSame(i).labelsMLIRM = labelsMLIRM(i);
    dataSame(i).classificationScoreIRM = classificationScoreIRM(i);
end

% SAVE THE CLASSIFIER MANUALLY IF ITS GOOD
eventClassifier_filespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/eventClassifier_RQSWDs.mat');
save(eventClassifier_filespec, 'trainedClassifier')

eventClassifier_filespec =  strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/eventClassifier_RQSWDs_trainingMat.mat');
save(eventClassifier_filespec, 'fullMat')

% PLOTS PLOTS PLOTS PLOTS PLOTS
predictorNames = fullMat.Properties.VariableNames(variablesvariables);

%{
% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsS1, 'br',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by S1');
print('~/Desktop/predictorVariables_S1.png', '-dpng');
close all

% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsS2, 'br',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by S2');
print('~/Desktop/predictorVariables_S2.png', '-dpng');
close all

% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsS3, 'rb',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by S3');
print('~/Desktop/predictorVariables_S3.png', '-dpng');
close all

% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsS4, 'rb',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by S4');
print('~/Desktop/predictorVariables_S4.png', '-dpng');
close all


% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsML, 'rb',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by SVM');
print('~/Desktop/predictorVariables_SVM.png', '-dpng');
close all

% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsHC, 'grb',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by Human Consensus');
print('~/Desktop/predictorVariables_ScorerConcensus.png', '-dpng');
close all

% plot each of the predictor variables against one another
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], labelsCL, 'rby',[],[],'on','stairs', predictorNames)
title('Predicotor Variables labeled by Clustering');
print('~/Desktop/predictorVariables_Clustering.png', '-dpng');
close all
%}

% INTER-RATER Agreement CAN IMPROVE THIS PLOT TO LOOK LIKE SOME OF THE ONES ABOVE
%{
figure('units', 'inch', 'pos', [10 10 20 20]) 
gplotmatrix(table2array(fullMat(:,variablesvariables)), [], weightedResponseWeights, 'ryg',[],[],'on','stairs', predictorNames);
title('Inter-Rater Agreement');
print('~/Desktop/predictorVariables_interraterAgreement.png', '-dpng');
close all
%}

%
% PLOT CONSISTENCY INTRA-RATER RELIABLITY
figure('units', 'inch', 'pos', [3 10 25 25]);
suptitle('Intra-rater Reliablity')
j = 1;
for i = 1:50

    theOnes = events(i,:);
    S1average(i) = mean(samelabelsS1(theOnes));
    S2average(i) = mean(samelabelsS2(theOnes));
    S3average(i) = mean(samelabelsS3(theOnes));
    S4average(i) = mean(samelabelsS4(theOnes));
    temp = round(S1average(i));
    propCorrectS1(i) = length(find(samelabelsS1(theOnes) == temp)) / nReps;
    temp = round(S2average(i));
    propCorrectS2(i) = length(find(samelabelsS2(theOnes) == temp)) / nReps;
    temp = round(S3average(i));
    propCorrectS3(i) = length(find(samelabelsS3(theOnes) == temp)) / nReps;
    temp = round(S4average(i));
    propCorrectS4(i) = length(find(samelabelsS4(theOnes) == temp)) / nReps;
    
    temp = dataSame(theOnes(1)).signalClips;
    temp = abs(cwt(temp, 256, 'amor'));
    temp = temp(:,1:lengthClip);
    if j < 50
        subplot(7, 7, j)
            clims = [0 3];
            imagesc(temp, clims)
            hold on;
            plot(1:lengthClip, (dataSame(theOnes(1)).signalClips(1:lengthClip) * -3) + 70, 'y');
            set(gca,'fontsize',15)
            xlim([0+512, 2048-512]);
            xlabel('Time (s)')
            axis off

            % add a patch for the windows that viewers saw when scoring seizures
            eventStart = 1024 - ceil(.5*dataSame(i).seizureDuration); % these need to be filled in.
            eventStop = 1024 + ceil(.5*dataSame(i).seizureDuration); % these need to be filled in.

            if round(S1average(i)) < 2
                displayS1 = 'notSWD';
                displayS1color = 'r';
            else
                displayS1 = 'SWD';  
                displayS1color = 'c';
            end

            if round(S2average(i)) < 2
                displayS2 = 'notSWD';
                displayS2color = 'r';
            else
                displayS2 = 'SWD';
                displayS2color = 'c';
            end

            if round(S3average(i)) < 2
                displayS3 = 'notSWD';
                displayS3color = 'r';
            else
                displayS3 = 'SWD';
                displayS3color = 'c';
            end

            if round(S4average(i)) < 2
                displayS4 = 'notSWD';
                displayS4color = 'r';
            else
                displayS4 = 'SWD';
                displayS4color = 'c';
            end

            if dataSame(theOnes(1)).labelsMLIRM == 2
                displaySVM = 'SWD';
                displaySVMcolor = 'c';
            else
                displaySVM = 'notSWD';
                displaySVMcolor = 'r';
            end    
            fontSize = 7;
            startLoc = 535;
            text(startLoc, 5, strcat({'S1: '}, num2str(propCorrectS1(i))), 'color', displayS1color, 'fontsize', fontSize)
            text(startLoc, 10, strcat({'S2: '}, num2str(propCorrectS2(i))), 'color', displayS2color, 'fontsize', fontSize)
            text(startLoc, 15, strcat({'S3: '}, num2str(propCorrectS3(i))), 'color', displayS3color, 'fontsize', fontSize)
            text(startLoc, 20, strcat({'S4: '}, num2str(propCorrectS4(i))), 'color', displayS4color, 'fontsize', fontSize)
            text(startLoc, 25, strcat({'SVM: '}, num2str(round(abs(dataSame(theOnes(1)).classificationScoreIRM), 2))), 'color', displaySVMcolor, 'fontsize', fontSize)
    end
    j = j + 1;

end

allSaverage = (S1average + S2average + S3average + S4average) ./ 4;
allPropCorrect = (propCorrectS1 + propCorrectS2 + propCorrectS3 + propCorrectS4) ./ 4;
print('~/Desktop/intraRater_Reliability_characterizeSWDs.png', '-dpng');
%print('~/Desktop/intraRater_Reliability_characterizeSWDs.pdf', '-dpdf', '-painters');

close all

% Average Interrater Agreement

figure('units', 'inch', 'pos', [10 10 15 15]) 
scatter3(fullMat.harmonic_std, fullMat.sixHz_max, fullMat.harmonic_mean, abs(classificationScore + .0001) * 10 , labelsMLColors)
xlabel('Std of the 16-32 Hz CWT Signal')
ylabel('Max of the 6Hz Hz CWT Signal')
zlabel('meanHarmonicSignal') 
set(gca,'fontsize',15)
title('Average Interrater Agreement')
hold on;
scatter3(interRaterMat.harmonic_std(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.sixHz_max(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.harmonic_mean(eventsMat(unique(eventsMat(1:50,1)),1)), 200 , allSaverage, 'filled')    
scatter3(interRaterMat.harmonic_std(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.sixHz_max(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.harmonic_mean(eventsMat(unique(eventsMat(1:50,1)),1)), 200 , 'k')     
colormap(flipud(redblue))
colorbar('Ticks',[1,1.5,2],'TickLabels',{'notSWD','Unsure','SWD'})
print('~/Desktop/averageInterraterAgreement.png', '-dpng');
close all

% Average Intrarater Reliability

figure('units', 'inch', 'pos', [10 10 15 15]) 
scatter3(fullMat.harmonic_std, fullMat.sixHz_max, fullMat.harmonic_mean, abs(classificationScore + .0001) * 10 , labelsMLColors)
xlabel('Std of the 16-32 Hz CWT Signal')
ylabel('Max of the 6Hz Hz CWT Signal')
zlabel('meanHarmonicSignal') 
title('Average Intrarater Reliability')
set(gca,'fontsize',15)
hold on;
scatter3(interRaterMat.harmonic_std(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.sixHz_max(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.harmonic_mean(eventsMat(unique(eventsMat(1:50,1)),1)), 200 , allPropCorrect, 'filled')    
scatter3(interRaterMat.harmonic_std(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.sixHz_max(eventsMat(unique(eventsMat(1:50,1)),1)), interRaterMat.harmonic_mean(eventsMat(unique(eventsMat(1:50,1)),1)), 200 , 'k')     
colormap(flipud(autumn))
colorbar()
print('~/Desktop/averageIntraraterReliability.png', '-dpng');
close all
%{

%% optimize parameters for linear kernel?

boxConstraint = [1, 2, 3, 4];
costValue = [1, 1.5, 2];

counter = 1;
for ii = costValue
    for jj = boxConstraint
        [trainedClassifier, validationAccuracy_mean(counter), validationAccuracy_std(counter)] = trainSWDClassifier_linear(humanMat, ii, jj);

        display(strcat({'Valdiation Accuracy with Training Data: '}, num2str(validationAccuracy_mean), {'+/-'}, num2str(validationAccuracy_std))); 

        % now predict on the full data set
        [labelsML, classificationScore] = trainedClassifier.predictFcn(fullMat);
        classificationScore = classificationScore(:,1); % the second colum is just the inverse of the first column when you use a gaussian kernel CHANGE IF CHANGING KERNEL


        labelsML = labelsML';

        % make color labels for labelsML and establish cal
        for i = 1:length(labelsML)
            if labelsML(i) == 2
                labelsMLColors(i,:) = [0, 0, 255]; % cyan rgb
            else
                labelsMLColors(i,:) = [255, 0, 0]; % red rgb

            end
        end

        % add this back to the data object for the inter-rater stuff later
        for i = 1:size(data, 2)
            data(i).labelsML = labelsML(i);
            data(i).classificationScore = classificationScore(i);
        end

        % assess quality of fit... and in relation to other scoreres
        labelslist = {'S1', 'S2', 'S3', 'S4', 'ML', 'HC', 'CC'};

        S1_S2 = sum(labelsS1 == labelsS2) / fullLength;
        S1_S3 = sum(labelsS1 == labelsS3) / fullLength;
        S1_S4 = sum(labelsS1 == labelsS4) / fullLength;
        S1_ML = sum(labelsS1 == labelsML) / fullLength;
        S1_HC = sum(labelsS1(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
        S1_CC = 1;
        S2_S3 = sum(labelsS2 == labelsS3) / fullLength;
        S2_S4 = sum(labelsS2 == labelsS4) / fullLength;
        S2_ML = sum(labelsS2 == labelsML) / fullLength;
        S2_HC = sum(labelsS2(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
        S2_CC = 1;
        S3_S4 = sum(labelsS3 == labelsS4) / fullLength;
        S3_ML = sum(labelsS3 == labelsML) / fullLength;
        S3_HC = sum(labelsS3(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
        S3_CC = 1;
        S4_ML = sum(labelsS4 == labelsML) / fullLength;
        S4_HC = sum(labelsS4(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
        S4_CC = 1;
        ML_HC(counter) = sum(labelsML(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
        ML_CC = 1;
        HC_CC = 1;

        matrixmatrix = [1,  S1_S2,  S1_S3,  S1_S4,  S1_ML,  S1_HC,  S1_CC; ...
                        0,  1,      S2_S3,  S2_S4,  S2_ML,  S2_HC,  S2_CC; ...
                        0,  0,      1,      S3_S4,  S3_ML,  S3_HC,  S3_CC; ...
                        0,  0,      0,      1,      S4_ML,  S4_HC,  S4_CC; ...
                        0,  0,      0,      0,      1,      ML_HC(counter),  ML_CC; ...
                        0,  0,      0,      0,      0,      1,      HC_CC;
                        0,  0,      0,      0,      0,      0,      1];

        sTable = array2table(matrixmatrix,'RowNames',labelslist,'VariableNames',labelslist) 

        % CALCULATE SOME RATIOS AND SHIT HERE
        averageHuman = mean([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
        stdHuman = std([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
        display(strcat({'Average Agreement Between Humans: '}, num2str(averageHuman), {'+/-'}, num2str(stdHuman)))
        averageML_Human(counter) = mean([S1_ML, S2_ML, S3_ML, S4_ML]);
        stdML_Human(counter) = std([S1_ML, S2_ML, S3_ML, S4_ML]);
        display(strcat({'Average Agreement Between ML and Humans: '}, num2str(averageML_Human(counter)), {'+/-'}, num2str(stdML_Human(counter))))


        % this wont work unless you build the dataSame structure down below already
        [labelsMLIRM classificationScoreIRM] = trainedClassifier.predictFcn(interRaterMat);
        classificationScoreIRM = classificationScoreIRM(:,1);


        lengthClip = 2048;
        buffer = 10;

        for i = 1:size(dataSame, 2)
            samelabelsS1(i) = dataSame(i).labelsS1;
            samelabelsS2(i) = dataSame(i).labelsS2;
            samelabelsS3(i) = dataSame(i).labelsS3;
            samelabelsS4(i) = dataSame(i).labelsS4;
            dataSame(i).labelsMLIRM = labelsMLIRM(i);
            dataSame(i).classificationScoreIRM = classificationScoreIRM(i);
        end

        for i = 1:length(dataSame)
            blahblahblah(i) = dataSame(i).labelsMLIRM;
            normScores2(i) = dataSame(i).classificationScoreIRM;
        end

        mdl = fitlm(normScores2(1:50),  1- (allSaverage - 1));

        RsquaredParam(counter) = mdl.Rsquared.Ordinary;
        mattsFuckingRatio(counter) = averageML_Human(counter) / averageHuman;
        costValueList(counter) = ii;
        boxConstraintList(counter) = jj;


        counter = counter + 1;
    end
end

costValueList = costValueList';
boxConstraintList = boxConstraintList';
averageML_Human = averageML_Human';
stdML_Human = stdML_Human';
ML_HC = ML_HC';
validationAccuracy_mean = validationAccuracy_mean';
validationAccuracy_std = validationAccuracy_std';
mattsFuckingRatio = mattsFuckingRatio';
RsquaredParam = RsquaredParam';

% print out the useful variables
optimizationParamsLinear = table(costValueList, boxConstraintList, averageML_Human, stdML_Human, ML_HC, validationAccuracy_mean, validationAccuracy_std, RsquaredParam);

parameters_filespec =  strcat(labDataDrive, 'jonesLab_data/sleep_and_seizures/EEG_data/RQ/optimizationParamsLinear.mat');
save(parameters_filespec, 'optimizationParamsLinear')

%% plot these parameters to make them useful
parameters_filespec =  strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/optimizationParamsLinear.mat');
load(parameters_filespec)

subplot(2,3, 1)
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3, '*')
    hold on;
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('SD of the 5x Cross Validation')
    colormap(redgreenblue) % XXX


subplot(2,3, 2)
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3, '*')
    hold on;
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('SD of the Agreement between the SVM and Human Scorers')
    colormap(redgreenblue) % XXX
    legend('no cost', 'cost = 1.5', 'cost = 2', 'cost = 2.5', 'cost = 3', 'Location', 'northwest')


    
subplot(2,3, 3)
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3, '*')
    hold on;
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.RsquaredParam(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('Agreement between the HC and SVM')
    colormap(redgreenblue) % XXX
    %colorbar
    
    
subplot(2,3, 4)
    scatter(optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3, '*')
    hold on;
    scatter(optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('SD of the 5x Cross Validation')
    ylabel('SD of the Agreement between the SVM and Human Scorers')


subplot(2,3, 5)
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3 ,'*')
    hold on;
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.stdML_Human(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('Agreement between the HC and SVM')
    ylabel('SD of the Agreement between the SVM and Human Scorers')

    
subplot(2,3, 6)
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1) * 3, '*')
    hold on;
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 1.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 1.5) * 3, '+')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2) * 3, 'o')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 2.5), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 2.5) * 3, 'v')
    scatter(optimizationParamsLinear.ML_HC(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.validationAccuracy_std(optimizationParamsLinear.costValueList == 3), optimizationParamsLinear.boxConstraintList(optimizationParamsLinear.costValueList == 3) * 3, '^')

    grid on;
    xlabel('Agreement between the HC and SVM')
    ylabel('SD of the 5x Cross Validation')
    %colorbar

print(gcf, '~/Desktop/Figure_optimization2Linear.pdf', '-dpdf', '-r300')
%}

%% optimize parameters for gaussian kernel
% set up a grid search for 3 parameters (box constraint, cost function, and kernel scale

kernelScale = [1, 5, 10, 15, 20];
boxConstraint = [5, 10, 15, 20, 25];
costValue = [1, 1.5, 2, 2.5, 3];

counter = 1;
for ii = costValue
    for jj = boxConstraint
        for kk = kernelScale
            [trainedClassifier, validationAccuracy_mean(counter), validationAccuracy_std(counter)] = trainSWDClassifier_gaussian(humanMat, ii, jj, kk);
 
            display(strcat({'Valdiation Accuracy with Training Data: '}, num2str(validationAccuracy_mean), {'+/-'}, num2str(validationAccuracy_std))); 

            % now predict on the full data set
            [labelsML, classificationScore] = trainedClassifier.predictFcn(fullMat);
            classificationScore = classificationScore(:,1); % the second colum is just the inverse of the first column when you use a gaussian kernel CHANGE IF CHANGING KERNEL


            labelsML = labelsML';

            % make color labels for labelsML and establish cal
            for i = 1:length(labelsML)
                if labelsML(i) == 2
                    labelsMLColors(i,:) = [0, 0, 255]; % cyan rgb
                else
                    labelsMLColors(i,:) = [255, 0, 0]; % red rgb

                end
            end

            % add this back to the data object for the inter-rater stuff later
            for i = 1:size(data, 2)
                data(i).labelsML = labelsML(i);
                data(i).classificationScore = classificationScore(i);
            end

            % assess quality of fit... and in relation to other scoreres
            labelslist = {'S1', 'S2', 'S3', 'S4', 'ML', 'HC', 'CC'};

            S1_S2 = sum(labelsS1 == labelsS2) / fullLength;
            S1_S3 = sum(labelsS1 == labelsS3) / fullLength;
            S1_S4 = sum(labelsS1 == labelsS4) / fullLength;
            S1_ML = sum(labelsS1 == labelsML) / fullLength;
            S1_HC = sum(labelsS1(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
            S1_CC = 1;
            S2_S3 = sum(labelsS2 == labelsS3) / fullLength;
            S2_S4 = sum(labelsS2 == labelsS4) / fullLength;
            S2_ML = sum(labelsS2 == labelsML) / fullLength;
            S2_HC = sum(labelsS2(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
            S2_CC = 1;
            S3_S4 = sum(labelsS3 == labelsS4) / fullLength;
            S3_ML = sum(labelsS3 == labelsML) / fullLength;
            S3_HC = sum(labelsS3(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
            S3_CC = 1;
            S4_ML = sum(labelsS4 == labelsML) / fullLength;
            S4_HC = sum(labelsS4(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
            S4_CC = 1;
            ML_HC(counter) = sum(labelsML(labelsHC == 1 | labelsHC == 2) == labelsHC(labelsHC == 1 | labelsHC == 2)) / sum(labelsHC == 1 | labelsHC == 2);
            ML_CC = 1;
            HC_CC = 1;

            matrixmatrix = [1,  S1_S2,  S1_S3,  S1_S4,  S1_ML,  S1_HC,  S1_CC; ...
                            0,  1,      S2_S3,  S2_S4,  S2_ML,  S2_HC,  S2_CC; ...
                            0,  0,      1,      S3_S4,  S3_ML,  S3_HC,  S3_CC; ...
                            0,  0,      0,      1,      S4_ML,  S4_HC,  S4_CC; ...
                            0,  0,      0,      0,      1,      ML_HC(counter),  ML_CC; ...
                            0,  0,      0,      0,      0,      1,      HC_CC;
                            0,  0,      0,      0,      0,      0,      1];

            sTable = array2table(matrixmatrix,'RowNames',labelslist,'VariableNames',labelslist) 

            % CALCULATE SOME RATIOS AND SHIT HERE
            averageHuman = mean([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
            stdHuman = std([S1_S2, S1_S3, S1_S4, S2_S3, S2_S4, S3_S4]);
            display(strcat({'Average Agreement Between Humans: '}, num2str(averageHuman), {'+/-'}, num2str(stdHuman)))
            averageML_Human(counter) = mean([S1_ML, S2_ML, S3_ML, S4_ML]);
            stdML_Human(counter) = std([S1_ML, S2_ML, S3_ML, S4_ML]);
            display(strcat({'Average Agreement Between ML and Humans: '}, num2str(averageML_Human(counter)), {'+/-'}, num2str(stdML_Human(counter))))


            % this wont work unless you build the dataSame structure down below already
            [labelsMLIRM classificationScoreIRM] = trainedClassifier.predictFcn(interRaterMat);
            classificationScoreIRM = classificationScoreIRM(:,1);


            lengthClip = 2048;
            buffer = 10;

            for i = 1:size(dataSame, 2)
                samelabelsS1(i) = dataSame(i).labelsS1;
                samelabelsS2(i) = dataSame(i).labelsS2;
                samelabelsS3(i) = dataSame(i).labelsS3;
                samelabelsS4(i) = dataSame(i).labelsS4;
                dataSame(i).labelsMLIRM = labelsMLIRM(i);
                dataSame(i).classificationScoreIRM = classificationScoreIRM(i);
            end

            for i = 1:length(dataSame)
                blahblahblah(i) = dataSame(i).labelsMLIRM;
                normScores2(i) = dataSame(i).classificationScoreIRM;
            end

            mdl = fitlm(normScores2(1:50),  1- (allSaverage - 1));

            RsquaredParam(counter) = mdl.Rsquared.Ordinary;
            mattsFuckingRatio(counter) = averageML_Human(counter) / averageHuman;
            costValueList(counter) = ii;
            boxConstraintList(counter) = jj;
            kernelScaleList(counter) = kk;


            counter = counter + 1;

        end
    end
end

costValueList = costValueList';
boxConstraintList = boxConstraintList';
kernelScaleList = kernelScaleList';
averageML_Human = averageML_Human';
stdML_Human = stdML_Human';
ML_HC = ML_HC';
validationAccuracy_mean = validationAccuracy_mean';
validationAccuracy_std = validationAccuracy_std';
mattsFuckingRatio = mattsFuckingRatio';
RsquaredParam = RsquaredParam';

% print out the useful variables
optimizationParams = table(costValueList, boxConstraintList, kernelScaleList, averageML_Human, stdML_Human, ML_HC, validationAccuracy_mean, validationAccuracy_std, RsquaredParam);

parameters_filespec =  strcat(labDataDrive, 'jonesLab_data/sleep_and_seizures/EEG_data/RQ/optimizationParams.mat');
save(parameters_filespec, 'optimizationParams')

% plot out some of the hyperparameters -- probably figure for paper.

% plot these parameters to make them useful
parameters_filespec =  strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/optimizationParams.mat');
load(parameters_filespec)

subplot(2,3, 1)
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1.5), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2.5), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 3), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('SD of the 5x Cross Validation')
    colormap(redgreenblue) % XXX


subplot(2,3, 2)
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1), optimizationParams.stdML_Human(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2), optimizationParams.stdML_Human(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 3), optimizationParams.stdML_Human(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('SD of the Agreement between the SVM and Human Scorers')
    colormap(redgreenblue) % XXX
    legend('no cost', 'cost = 1.5', 'cost = 2', 'cost = 2.5', 'cost = 3', 'Location', 'northwest')


    
subplot(2,3, 3)
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1), optimizationParams.ML_HC(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 1.5), optimizationParams.ML_HC(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2), optimizationParams.ML_HC(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 2.5), optimizationParams.ML_HC(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.RsquaredParam(optimizationParams.costValueList == 3), optimizationParams.ML_HC(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('R^2 between SVM and Repeat Events')
    ylabel('Agreement between the HC and SVM')
    colormap(redgreenblue) % XXX
    %colorbar
    
    
subplot(2,3, 4)
    scatter(optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1), optimizationParams.stdML_Human(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2), optimizationParams.stdML_Human(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 3), optimizationParams.stdML_Human(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('SD of the 5x Cross Validation')
    ylabel('SD of the Agreement between the SVM and Human Scorers')


subplot(2,3, 5)
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 1), optimizationParams.stdML_Human(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 1.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 2), optimizationParams.stdML_Human(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 2.5), optimizationParams.stdML_Human(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 3), optimizationParams.stdML_Human(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('Agreement between the HC and SVM')
    ylabel('SD of the Agreement between the SVM and Human Scorers')

    
subplot(2,3, 6)
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 1), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1), '*')
    hold on;
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 1.5), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 1.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 1.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 1.5), '+')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 2), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2), 'o')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 2.5), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 2.5), optimizationParams.boxConstraintList(optimizationParams.costValueList == 2.5) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 2.5), 'v')
    scatter(optimizationParams.ML_HC(optimizationParams.costValueList == 3), optimizationParams.validationAccuracy_std(optimizationParams.costValueList == 3), optimizationParams.boxConstraintList(optimizationParams.costValueList == 3) * 3, optimizationParams.kernelScaleList(optimizationParams.costValueList == 3), '^')

    grid on;
    xlabel('Agreement between the HC and SVM')
    ylabel('SD of the 5x Cross Validation')
    %colorbar

print(gcf, '~/Desktop/Figure_optimization2.pdf', '-dpdf', '-r300')

%
%% plots for the manuscript

% variable for the weighting of yes/no events for the original scoring.. only includes unique events for display purposes (not for anlyses) 
uniqueInd = zeros(length(labelsS1), 1);
[a, b] = unique(fullMat);
uniqueInd(b) = 1;
keepInd = labelsS1 ~=0 & labelsS2 ~=0 & labelsS3 ~=0 & labelsS4 ~=0 & uniqueInd' ~= 0;
newWeights = mean([labelsS1(keepInd)', labelsS2(keepInd)', labelsS3(keepInd)', labelsS4(keepInd)']');
newWeights = newWeights - 1;

altColors = flipud(colormap(autumn));
colors = colormap(redgreenblue());
%colors = colormap(flipud(cbrewer('div', 'Spectral', 64)));

% variable name for plotting
publicationVariables = {'4.4-8.2 Hz Mean', '4.4-8.2 Hz Std', 'Max CWT Value at ~ 4.4-8.2 Hz', '8.8-16.4 Hz Mean', '8.8-16.4 Hz Std', 'Max CWT Value at ~ 8.8-16.4 Hz', '17.6-32.8 Hz Mean', '17.6-32.8 Hz Std', 'Max CWT Value at ~ 17.6-32.8 Hz', '35.1-65.5 Hz Mean', '35.1-65.5 Hz Std', 'Max CWT Value at ~ 35.1-65.5 Hz'};
%%
% ----- Fig: 1 CWT image
figure('units', 'inch', 'pos', [10 10 12 3]) 
plot(data(291).signalClips, 'k'); % 236
xlim([1,2048])
ylim([-8, 8])
set(gca, 'fontsize', 15)
[sx, sy] = scalebars(gca, [1400, 5], [25, 31], {'Sample Points = 1 s', 'Norm. EEG Amp.'}, 'arial', 15, 2);
axis off
set(gcf, 'PaperSize', [20 20]) 
print(gcf, '~/Desktop/Figure_1a.pdf', '-dpdf', '-r300')
figure('units', 'inch', 'pos', [10 10 13 4]) 
cwt(data(291).signalClips, 'amor', 256)
colormap(colors)

set(gca, 'fontsize', 15)
set(gcf, 'PaperSize', [20 20]) 
print(gcf, '~/Desktop/Figure_1b.pdf', '-dpdf', '-r300')

%% New Fig 2
figure('units', 'inch', 'pos', [10 10 20 20]) 
temp = [1,2,3,7,8,9];
gplotmatrix(table2array(fullMat(:,temp)), [], labelsHC, 'gbr',[],[],'on','stairs', predictorNames(temp))
title('Predicotor Variables labeled by Human Consensus');

    set(gcf, 'PaperSize', [20 20]) 
    print(gcf, '~/Desktop/newFigure_2.pdf', '-painters', '-dpdf', '-r300')

%%
% ----- Fig 2: 8 example events and these events along with the 2050 events that were scored by reviers
figure('units', 'inch', 'pos', [10 10 15 15]) 
yLimits = [-8, 8];

%colormap(altColors)
colormap(colors)


% example events to display in this figure
eventNumbers = [9, 20, 3, 6, 44, 10, 8, 50];

subplot(6, 12, 1:3)
    eventNum = eventNumbers(1);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '1', 'fontsize', 20)

    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end

    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 4:6)
    eventNum = eventNumbers(2);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '2', 'fontsize', 20)
    
    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end

    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
     
subplot(6, 12, 7:9)
    eventNum = eventNumbers(3);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
     text(600, 6, '3', 'fontsize', 20)
   % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 10:12)
    eventNum = eventNumbers(4);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '4', 'fontsize', 20)
    [sx, sy] = scalebars(gca, [1400, 5], [25, 31], {'Sample Points = 1 s', 'Norm. EEG Amp.'}, 'arial', 15, 2);

    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 13:15)
    eventNum = eventNumbers(5);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '5', 'fontsize', 20)
    
    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 16:18)
    eventNum = eventNumbers(6);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '6', 'fontsize', 20)
    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 19:21)
    eventNum = eventNumbers(7);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '7', 'fontsize', 20)
    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
subplot(6, 12, 22:24)
    eventNum = eventNumbers(8);
    theOnes = events(eventNum,:);
    plot(1:lengthClip, dataSame(theOnes(1)).signalClips(1:lengthClip), 'k')
    xlim([512, 1536]);
    ylim(yLimits)
    axis off;
    text(600, 6, '8', 'fontsize', 20)

    % add a patch for the windows that viewers saw when scoring seizures
    colorNum = (allSaverage(eventNum) - 1);
    if colorNum == 0
        colorcolor = colors( 1,:);
    else
        colorcolor = colors( floor(64 * colorNum) ,:);
    end
 
    eventStart = 1024 - ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    eventStop = 1024 + ceil(.5*dataSame(eventNum).seizureDuration); % these need to be filled in.
    patch([eventStart - buffer, eventStop + buffer, eventStop + buffer, eventStart - buffer], [-7, -7, -10, -10], colorcolor, 'FaceColor', colorcolor, 'EdgeColor', 'none')
    
% setups for the new few plots
variableA = 3;
variableB = 6;
variableC = 9;
variableD = 12;

% intrarater agreement all events  

pt = .40;

subplot(6, 12, [25:28, 37:40])
    scatter(table2array(fullMat((keepInd),variableA)), table2array(fullMat((keepInd),variableC)), 10, newWeights, 'filled', 'MarkerFaceAlpha',pt,'MarkerEdgeAlpha',pt)
    hold on;
        set(gca, 'fontsize', 15)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(1),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(1),1)),1)), '1', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(2),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(2),1)),1)), '2', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(3),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(3),1)),1)), '3', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(4),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(4),1)),1)), '4', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(5),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(5),1)),1)), '5', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(6),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(6),1)),1)), '6', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(7),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(7),1)),1)), '7', 'fontsize', 20)
    text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(8),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(8),1)),1)), '8', 'fontsize', 20)

% interrater agreement (50 events)
subplot(6, 12, [29:32, 41:44])
    scatter(table2array(fullMat((keepInd),variableA)), table2array(fullMat((keepInd),variableC)), 10, newWeights, 'filled', 'MarkerFaceAlpha',pt,'MarkerEdgeAlpha',pt)
    hold on;
    set(gca, 'fontsize', 15)
    scatter(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(:,1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(:,1)),1)), 100 , allSaverage - 1, 'filled')    
    scatter(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(:,1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(:,1)),1)), 100 , 'k')  
    
% intrarater agreement (50 events)
subplot(6, 12, [33:36, 45:48])
    scatter(table2array(fullMat((keepInd),variableA)), table2array(fullMat((keepInd),variableC)), 10, newWeights, 'filled', 'MarkerFaceAlpha',pt,'MarkerEdgeAlpha',pt)
    hold on;
    set(gca, 'fontsize', 15)
    scatter(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(:,1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(:,1)),1)), 100 , (allPropCorrect - min(allPropCorrect)) / ( max(allPropCorrect) - min(allPropCorrect) ), 'filled') % norm_data = (allPropCorrect - min(allPropCorrect)) / ( max(allPropCorrect) - min(allPropCorrect) )
    scatter(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(:,1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(:,1)),1)), 100 , 'k')   
    colorbar()

% scorer 1 labels
subplot(6, 12, [49:52, 61:64])

    scatter(table2array(fullMat(labelsS1 ==1 ,variableA)), table2array(fullMat(labelsS1 == 1,variableC)), 10, colors(1,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    hold on
    set(gca, 'fontsize', 15)
    scatter(table2array(fullMat(labelsS1 ==2,variableA)), table2array(fullMat(labelsS1 == 2,variableC)), 10, colors(end,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    ylabel(publicationVariables{variableC})

% scorer 2 labels
subplot(6, 12, [53:56, 65:68])
    scatter(table2array(fullMat(labelsS2 ==1,variableA)), table2array(fullMat(labelsS2 == 1,variableC)), 10, colors(1,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    hold on
    set(gca, 'fontsize', 15)
    scatter(table2array(fullMat(labelsS2 ==2,variableA)), table2array(fullMat(labelsS2 == 2,variableC)), 10, colors(end,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    xlabel(publicationVariables{variableA})

% scorer 3
subplot(6, 12, [57:60, 69:72])
    scatter(table2array(fullMat(labelsS3 ==1,variableA)), table2array(fullMat(labelsS3 == 1,variableC)), 10, colors(1,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    set(gca, 'fontsize', 15)
    hold on
    scatter(table2array(fullMat(labelsS3 ==2,variableA)), table2array(fullMat(labelsS3 == 2,variableC)), 10, colors(end,:), 'filled', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25)
    set(gcf, 'PaperSize', [20 20]) 
    print(gcf, '~/Desktop/Figure_2.pdf', '-dpdf', '-r300')

    %%
% ----- Fig 3:

% normalize the scores from the svm classification so everything is between 0 and 1.
normScores= (classificationScore - min(classificationScore)) / (max(classificationScore) - min(classificationScore));
normScores = normScores(keepInd);
figure('units', 'inch', 'pos', [10 10 10 5]) 
colormap(colors)

subplot(1, 8, 1:4)
    scatter(table2array(fullMat((keepInd),variableA)), table2array(fullMat((keepInd),variableC)), 10, 1 - normScores, 'filled', 'MarkerFaceAlpha',.25,'MarkerEdgeAlpha',.25)
    set(gca, 'fontsize', 15)
    hold on;
    for i = 1:length(eventNumbers)
        if dataSame(eventNumbers(i)).labelsMLIRM == 2
            eventColor(i,:) = colors(end,:);
        else
            eventColor(i,:) = colors(1,:);
        end
    end
    num = 1; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 2; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 3; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 4; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 5; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 6; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 7; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    num = 8; text(interRaterMat.sixHz_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), interRaterMat.harmonic_max(eventsMat(unique(eventsMat(eventNumbers(num),1)),1)), num2str(num), 'fontsize', 20, 'color', eventColor(num,:))
    xlabel(publicationVariables{variableA})
    %colorbar()


subplot(1, 8, 5:8)
    for i = 1:length(dataSame)
        blahblahblah(i) = dataSame(i).labelsMLIRM;
        normScores2(i) = dataSame(i).classificationScoreIRM;
    end
    noInd = find(blahblahblah(1:50) ==  1);
    yesInd = find(blahblahblah(1:50) ==  2);
    mdl = fitlm(normScores2(1:50),  1- (allSaverage - 1));

    plot(mdl);
    set(gca, 'fontsize', 15)
    xlabel('Norm. Dist. to SV');
    ylabel('Avg. Rater Classification')
    legend off;
    ylim([0, 1])
    title('')
    text(.5, .1, strcat({'R^2: '}, num2str(round(mdl.Rsquared.Ordinary, 3)), {' P < 0.001'}), 'fontsize', 12)
    set(gcf, 'PaperSize', [20 20]) 
    print(gcf, '~/Desktop/Figure_3b.pdf', '-dpdf', '-r300')

%% ----- Fig 4: Agreement Matrix
sTable(1:7,1:7)
display(strcat({'Average Agreement Between Humans: '}, num2str(averageHuman), {'+/-'}, num2str(stdHuman)))
display(strcat({'Average Agreement Between ML and Humans: '}, num2str(averageML_Human), {'+/-'}, num2str(stdML_Human)))

% ----- Fig 5 found in summarizeIndividualRecords.m

% ----- Fig 6 found in summarizeIndividualRecords.m

%%

%{
% Regress PLOT INTRARATER RELIABLITY AGAINST IMPORTANT VARIABLES

figure('units', 'inch', 'pos', [10 10 12 12]) 

plot(abs(classificationScore(round(allSaverage) == 2)), jitter(allPropCorrect(round(allSaverage) == 2)), '.c', 'markersize', 25)
hold on;
plot(abs(classificationScore(round(allSaverage) == 1)), jitter(allPropCorrect(round(allSaverage) == 1)), '.r', 'markersize', 25)
%legend({'SWD', 'notSWD'}, 'Location','SouthEast')
ylabel('Intra-rater Reliability (Prop.)')
xlabel('SVM Classification Score')
set(gca, 'fontsize', 18)
ylim([0.741, 1.009])

[ryes, m, b] = regression(abs(classificationScore(round(allSaverage) == 2)) , allPropCorrect(round(allSaverage) == 2));
x = 0:0.1:round(max(abs(classificationScore)));
y = m*x + b;
line(x, y, 'color', 'b', 'linewidth', 2)

[rno, m, b] = regression(abs(classificationScore(round(allSaverage) == 1)) , allPropCorrect(round(allSaverage) == 1));
x = 0:0.1:round(max(abs(classificationScore)));
y = m*x + b;
line(x, y, 'color', 'r', 'linewidth', 2)

text(2, .782, strcat({'SWD, R^2 = '}, num2str(ryes)), 'fontsize', 20, 'color', 'b')
text(2, .77, strcat({'notSWD, R^2 = '}, num2str(rno)), 'fontsize', 20, 'color', 'r')

print('~/Desktop/intraRater_Reliability_vs_classification_characterizeSWDs.png', '-dpng');
close all

% same thing as above plot but with std. harmonic signal rather than classification Probability

figure('units', 'inch', 'pos', [10 10 12 12]) 
temp = eventsMat(unique(eventsMat(:,1)),1);
temp(round(allSaverage) == 2)
plot(interRaterMat.stdHarmonicSignal(temp(round(allSaverage) == 2)), jitter(allPropCorrect(round(allSaverage) == 2)), '.c', 'markersize', 25)
hold on;
plot(interRaterMat.stdHarmonicSignal(temp(round(allSaverage) == 1)), jitter(allPropCorrect(round(allSaverage) == 1)), '.r', 'markersize', 25)
%legend({'SWD', 'notSWD'})
ylabel('Intra-rater Reliability (Prop.)')
xlabel('Std. Harmonic Signal')
set(gca, 'fontsize', 18)
ylim([0.741, 1.009])

[ryes, m, b] = regression(interRaterMat.stdHarmonicSignal(temp(round(allSaverage) == 2))' , allPropCorrect(round(allSaverage) == 2));
x = 0:0.1:350;
y = m*x + b;
line(x, y, 'color', 'b', 'linewidth', 2)

[rno, m, b] = regression(interRaterMat.stdHarmonicSignal(temp(round(allSaverage) == 1))' , allPropCorrect(round(allSaverage) == 1));
x = 0:0.1:350;
y = m*x + b;
line(x, y, 'color', 'r', 'linewidth', 2)

text(10, .782, strcat({'SWD, R^2 = '}, num2str(ryes)), 'fontsize', 20, 'color', 'b')
text(10, .77, strcat({'notSWD, R^2 = '}, num2str(rno)), 'fontsize', 20, 'color', 'r')

print('~/Desktop/intraRater_Reliability_vs_classification_characterizeSWDs2.png', '-dpng');
close all

%}

%{
figure('units', 'inch', 'pos', [10 10 12 12]) 
temp = eventsMat(unique(eventsMat(:,1)),1);
temp(round(allSaverage) == 2)
plot(interRaterMat.eventLength(temp(round(allSaverage) == 2)), jitter(allPropCorrect(round(allSaverage) == 2)), '.c', 'markersize', 25)
hold on;
plot(interRaterMat.eventLength(temp(round(allSaverage) == 1)), jitter(allPropCorrect(round(allSaverage) == 1)), '.r', 'markersize', 25)
legend({'SWD', 'notSWD'})
ylabel('Intra-rater Reliability (Prop.)')
xlabel('Event Length')
set(gca, 'fontsize', 18)

[r, m, b] = regression(interRaterMat.eventLength(temp(round(allSaverage) == 2))' , allPropCorrect(round(allSaverage) == 2));
x = 0:0.1:350;
y = m*x + b;
line(x, y, 'color', 'b', 'linewidth', 2)

[r, m, b] = regression(interRaterMat.eventLength(temp(round(allSaverage) == 1))' , allPropCorrect(round(allSaverage) == 1));
x = 0:0.1:350;
y = m*x + b;
line(x, y, 'color', 'r', 'linewidth', 2)

%}

% generate interRater predictor variables for event classifier and then plot them to see if they fill the same space as the rest of the data.

%



