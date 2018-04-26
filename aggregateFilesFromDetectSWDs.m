%% ----- Step 2 Aggregate Data 

% ----- import the namefile
namefilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleep_and_seizures_RQ_namefile.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'

% set sample rate
Fs = 256; % samples per seconde
epochDuration = 4; % seconds
epochPoints = Fs * epochDuration;

% summarize SWD counts and find exact peaks for each animal.
for i = 1:length(F.itemparams)
    for l = 1:length(F.itemparams{1, i}.recordDay)
        j = str2num(F.itemparams{1, i}.recordDay{1, l});
        animalID = F.itemparams{1, i}.animalID{1};
        sleepScorer = F.itemparams{1, i}.sleepScore{1,l};
        %genotype = F.itemparams{1, i}.genotype{1};
        
        % import seizure scores
        seizureScoreSpec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/', animalID, '_', num2str(j), '_exactSeizureLocations.mat'); % local for now, and just day 4
        if exist(seizureScoreSpec, 'file')

        seizureScore = load(seizureScoreSpec);

            % restrict to seizures happening within a 24 hour window
            clear within24
            win24 = 22118400; % this is 24 hours in epoch points and the largest point in an EEG record that we should examine since there are overlapping points.
            if size(seizureScore.output, 2) > 0 
                for ii = 1:size(seizureScore.output, 2)
                    if seizureScore.output(1,ii).SWDLOCS(1) < win24
                        animal(i).day(j).SWDs(ii).SWDLOCS = seizureScore.output(1,ii).SWDLOCS;
                        animal(i).day(j).SWDs(ii).SWDPKS = seizureScore.output(1,ii).SWDPKS;
                        animal(i).day(j).seizureLength(ii) = range(seizureScore.output(1,ii).SWDLOCS) / Fs; % in seconds
                        animal(i).day(j).averageSeizureMag(ii) = mean(seizureScore.output(1,ii).SWDPKS); % in first differetial peak so not quite right. 
                        animal(i).day(j).stdSeizureMag(ii) = std(seizureScore.output(1,ii).SWDPKS);
                        animal(i).day(j).nSWDPKS(ii) = length(seizureScore.output(1,ii).SWDPKS);
                        animal(i).day(j).averagePeakWidth(ii) = mean(diff(seizureScore.output(1,ii).SWDLOCS)) / Fs; % in seconds
                        animal(i).day(j).stdPeakWidth(ii) = std(diff(seizureScore.output(1,ii).SWDLOCS)) / Fs; % in seconds
                        animal(i).day(j).classificationResponse(ii) = seizureScore.output(1,ii).classificationResponse;
                        animal(i).day(j).classificationScore(ii) = abs(seizureScore.output(1,ii).classificationScore);
                    end
                end
                predicoterMatrixfilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/', animalID, '_', num2str(j), '_predictorMatrix.mat');
                if exist(predicoterMatrixfilespec, 'file')
                    load(predicoterMatrixfilespec);
                    animal(i).day(j).predictorMatrix = predictorMat(1:size(seizureScore.output, 2),:);
                end
            end
        end
        % load the EDF file
        edffilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/EDFs/', animalID, '_', num2str(j), '.edf');
        animal(i).day(j).edffilespec = edffilespec;     
        % load sleep score file
        sleepscorefilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleepScoring/', sleepScorer, '/', animalID, '_', num2str(j), '.txt');
        if exist(sleepscorefilespec, 'file')
            animal(i).day(j).sleepScore = import_scores(sleepscorefilespec);
        end   
        

    end
end

% ----- save aggregated data

save(strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/aggregated_RQ_Data.mat'), 'animal');
