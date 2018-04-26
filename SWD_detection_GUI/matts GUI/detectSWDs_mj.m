function [ ] = detectSWDs( edffilespec )
%UNTITLED2 Summary of this function goes here
% This script identifies SWD epochs, presents them to a reviewer and allows
% the reviewer to identify the exact locations of the SWDs if they are
% present.
%
% Click on the first discernable peak of each seizure, then
% then the last peak. If you make a mistake then click the space bar to go
% backwards and repeat a record. Files are saved to the same folder 
%
% example:
%
% [basePath, ~, cookieMonster] = getUserPath();
% edffilespec = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/EDFs/Ronde_4.edf');
% detectSWDs( edffilespec);

% variable argument defaults

[basePath, ~, ~] = getUserPath();

% load EDF file
EDF = read_EDF_mj(edffilespec);


% calculate predictive matrix
disp('Generating the Matrix of Predictors...');
predictorMatrix = generateSWDPredictors( EDF );
% predictorMatrix = generateSWDPredictors_mj( EDF );

% load SVM classifier
load(strcat(basePath, '/jonesLab_git/jonesLab_code/Projects/sleep_and_seizures/EEG_Analysis/seizureIdentification/SVMClassifier_SWDs/polynomialSWDSVMClassifier.mat'));

% classify the SWDs
[yfit, ~] = trainedClassifier.predictFcn(predictorMatrix); 



%% seizure viewer

%key = 1;
loopStart = 1;
i = 1;
signal = EDF.D.edf.signalMat(:, 2);
sampleRate = 1024; % per 4 second epoch
recordCmp.edfLocs = find(yfit) * sampleRate;
figure('units', 'inch', 'pos', [1 5 20 5]);

  
while i < length(recordCmp.edfLocs)
    for i = loopStart:length(recordCmp.edfLocs)
        tagStart = recordCmp.edfLocs(i);
            if tagStart < 3 * sampleRate
                recordCmp.exactSeizureStart(i) = 0; % place start locations in object
                recordCmp.exactSeizureStop(i) = 0;
            else
                beforeWin = sampleRate * 2 - 1;
                afterWin = sampleRate * 2 - 1;
                seq = (tagStart - beforeWin):(tagStart + afterWin);
                plot(signal(seq));
                minimum = min(signal(seq)) + (min(signal(seq)) * 0.1); % minimum value for the EEG seq we are looking at plus 10 percent
                maximum = max(signal(seq)) + (max(signal(seq)) * 0.1); % same but for max
                bottom = -1 * max(abs([minimum, maximum])); % calculates where the bottom of the rectange that highlights the active epoch should go
                top = 2 * abs(bottom); % determines the top of the rectange
                rectangle('Position', [sampleRate, bottom, sampleRate, top], 'EdgeColor', 'r'); % places a rectangle on the plot that indicates the orignial tag location
                set(gca, 'xlim', [0, sampleRate*3]);
                
                % plot some labels and information
                text(1100, abs(bottom) - 0.01, ['record ', num2str(i), ' of ', num2str(length(recordCmp.edfLocs))]);
                %text(1100, bottom + .01, ['Distance from separator: ',
                %num2str(abs(distanceFromClassifierThresh(i, 1)))]); % not
                %sure why this is not represetitive of the distance from
                %the threshold, but its not..

                %waitforbuttonpress;
                [a, ~, key] = ginput(2); % code that allows user to input 2 clicks for start and stop
                recordCmp.exactSeizureStart(i) = round(a(1,1)) + (tagStart - beforeWin); % place start locations in object
                recordCmp.exactSeizureStop(i) = round(a(2,1)) + (tagStart - beforeWin); % place stop locations in object


                % go backwards if you press the spacebar 2x
                if key ~= 1
                    loopStart = i - 1;
                    break
                end
            end
    end
end

close all
  
%% export SWDFile

[a, b] = fileparts(edffilespec);
ds = cell2dataset(num2cell([recordCmp.exactSeizureStart; recordCmp.exactSeizureStop]'), 'VarNames', {'startSWD', 'endSWD'}); 
export(ds,'file', char(strcat(a, filesep, b, '_exactSeizureLocations.csv')),'delimiter',',')

end

