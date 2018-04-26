%% ----- Step 1 Generate Predictor Matricies and detectSWDs from SVM, and aggregate DATA

% ----- import the namefile
namefilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleep_and_seizures_RQ_namefile.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'

% ----- generate predictor matricies and detectSWDS

tic
for i = 1:length(F.itemparams)
    for l = 1:length(F.itemparams{1, i}.recordDay)
        j = str2num(F.itemparams{1, i}.recordDay{1, l});
        if str2num(F.itemparams{1, i}.SWDdetectOK{1, l}) == 1 
            animalID = F.itemparams{1, i}.animalID{1};
            edffilespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/EDFs/', animalID, '_', num2str(j), '.edf');
            eventClassifier_filespec = strcat(labDataDrive, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/eventClassifier_RQSWDs.mat');
            fs = 256;
            channel = 2;
            detectSWDs_automated(edffilespec, eventClassifier_filespec, fs, channel);
        end
    end
end
toc

% aggregate the data for summary output
aggregateFilesFromDetectSWDs

% summarize sleep and seizures in 24 hour displays
summarizeIndividualRecords

% summarize individual events for each plot
proofEvents