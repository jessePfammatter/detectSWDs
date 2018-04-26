% ----- add appropriate edf data (only snipits from SWD + display window
namefilespec = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleep_and_seizures_RQ_namefile.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'
load(strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/aggregated_RQ_Data.mat'));

% for each of the following animals_day combinations go download all the important clips of edf files
% animal list code 
animalLookup = [1, 4; 4, 4; 4, 5; 4, 6; 5, 4; 5, 5; 6, 4; 6, 5; 6, 6; 8, 2; 7, 1; 7, 2; 7, 3; 8, 3; 9, 2; 9, 3]; % all these animals I've done manual checks
clear data
zz = 1;

for i = 1:size(animalLookup, 1)
    EDF = read_EDF_mj(animal(animalLookup(i, 1)).day(animalLookup(i, 2)).edffilespec);
    signal = EDF.D.edf.signalMat(:,1:4);
    for ii = 1:length(animal(animalLookup(i, 1)).day(animalLookup(i, 2)).SWDs)
        
        
        seizureStart = animal(animalLookup(i, 1)).day(animalLookup(i, 2)).SWDs(ii).SWDLOCS(1);
        seizureStop = animal(animalLookup(i, 1)).day(animalLookup(i, 2)).SWDs(ii).SWDLOCS(end);
        center = floor(seizureStop - seizureStart);
        windowStart = floor(seizureStart + ceil(.5*center) - 1024);
        windowStop = floor(seizureStop - ceil(.5*center) + 1024);
        animal(animalLookup(i, 1)).day(animalLookup(i, 2)).SWDSignalClips{ii} = signal(windowStart:windowStop,:);
        
        % add to new data set
        data(zz).seizureStart = seizureStart;
        data(zz).seizureStop = seizureStop;
        data(zz).seizureDuration = center;
        data(zz).signalClips = signal(windowStart:windowStop,:);
        data(zz).animalInfo = animal(animalLookup(i, 1)).day(animalLookup(i, 2)).edffilespec;
        data(zz).originalSeizeNumber = ii;
        
        % increment counter
        zz = zz + 1;
        
        
    end
end

% reduce size of data set
dataOrig = data;
nStart = size(data, 2);
nSubsample = 2050;
subsampleInd = randsample(nStart, nSubsample);
data = data(subsampleInd);
dataOrig = dataOrig(subsampleInd);

%% calculate the variables used for the SVM so that you can select 50 events evenly across the n dimensional space.
predictors = generatePredictorsForEventClassifier(dataOrig);
%predictors = generatePredictorsForEventClassifier(D.data);

includedPredictorNames = predictors.Properties.VariableNames([false false true false true true false false false false false false true false true true false true false true]);
predictors = predictors(:,includedPredictorNames);

%% randomly select 50 swds to replicate 10x each
nDuplicate = 50;
duplicateList = randsample(nSubsample, nDuplicate);

%% check how evently the samples are distributed if you just select at random.

predictorsMat = table2array(predictors);

plot3(predictorsMat(duplicateList, 2), predictorsMat(duplicateList, 5), predictorsMat(duplicateList, 6), 'or');
hold on;
plot3(predictorsMat(:, 2), predictorsMat(:, 5), predictorsMat(:, 6), '.b');
%plot3(predictorsMat(:, 2), predictorsMat(:, 5), predictorsMat(:, 6), '.');
xlabel('Event Duration')
ylabel('Mean Harmonic Signal')
zlabel('Std Harmonic Signal')
grid on


%% if good then integrate them, shuffle the data, and build the data set for the SWDGUI

% a bunch from Ronde_4 and then a bunch of randoms
%duplicateList = [246, 247, 249, 250, 251, 254, 257, 287, 288, 292, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500];
zz = size(data, 2) + 1;
clear i
for j = 1:9
    for i = 1:length(duplicateList)
        ii = duplicateList(i);
        data(zz) = data(ii);
        zz = zz + 1;
    end
end

% shuffle the data

shuffledList = randperm(size(data, 2));
for i = 1:size(data, 2)
    data(i).presentationOrder = shuffledList(i);
end



%
D.data = data;
for i = 1:size(D.data, 2)
    D.data(i).responses = 0;
    D.data(i).comments = [];
end
D.data
D.zz = 1;


%% save the data

%
save('~/Desktop/SWDDATA_07:20:17.mat', 'D');
%load('/Users/jesse/Desktop/SWDDATA_07-20-17.mat');

