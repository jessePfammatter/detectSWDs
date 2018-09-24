% This file is designed to compare scoring between humans and the support
% vector machine (SVM) as characterized in the manuscript, 'An automated, machine learning-based detection algorithm for spike-wave discharges (SWDs) in a mouse model of absence epilepsy'

%% setup workspace

% this is used by the jonesLab to setup our path
[basePath, merlinPath, labDrive, user] = getUserPath();

seizurePath = '/gamma2R43Q/seizureScoring/';
edfPath = '/gamma2R43Q/EDFs/';
compareScoresOutputPath = '/gamma2R43Q/compareScores/';

% prepare table to compare records

% ronde_4 = rq, ronde_5 = ganx, sherlock 4 = RR, sherlock 6 = rr and ganx, tiki_4 = rq, tiki_5 = rq + ganx, vincent_4 = RQ, vincent_10 = rq and ganx)
files = cellstr({'Ronde_4'; 'Ronde_5'; 'Sherlock_6'; 'Tiki_4'; 'Vincent_4'});

% set some standards for the data (everything is in 256 Hz sample rate
epochDurationSampleRate = 1024;
epochDurationSeconds = 4;

%% start collecting the relevant data

fileStart = 1; % this is a counter to help build the appropriate files. The filepaths in the section below should be changed to wherever you have stored the relevant files
for file = fileStart:length(files)
  
  % filepath for each set of files
   temp{file}.edffilespec   = char(strcat(labDrive, edfPath, files(file), '.edf'));
   temp{file}.s1filespec    = char(strcat(labDrive, seizurePath, 'rm/', files(file), '.txt'));
   temp{file}.s2filespec 	= char(strcat(labDrive, seizurePath, 'jp/', files(file), '.tag'));
  
  % import edf file
  [temp{file}.edf.header, temp{file}.edf.signalHeader  temp{file}.edf.signalCell] = blockEdfLoad(temp{file}.edffilespec);
  
  % import all the scores
  temp{file}.s1 = import_scores(temp{file}.s1filespec);
  temp{file}.s2 = import_scores(temp{file}.s2filespec);
  
  % label seizure tags for each scorer
  temp{file}.s1label = 'SWD';
  temp{file}.s2label = 'SWD';
  
  % make a vector of zeros and ones to directly compare seizure scoring
  temp{file}.s1EventIndicator = strcmp(temp{file}.s1.score, temp{file}.s1label);
  temp{file}.s2EventIndicator = strcmp(temp{file}.s2.score, temp{file}.s2label);
  
  % index of
  temp{file}.s1EventIndex = find(temp{file}.s1EventIndicator);
  temp{file}.s2EventIndex = find(temp{file}.s2EventIndicator);

  % make a united list of all the epochs
  temp{file}.interEventIndex   = intersect(temp{file}.s1EventIndex, temp{file}.s2EventIndex);
  temp{file}.unionEventIndex   = union(temp{file}.s1EventIndex, temp{file}.s2EventIndex);
  
  recordCmp{file}           = array2table(temp{file}.s1.time(temp{file}.unionEventIndex), 'VariableNames', {'time'});
  recordCmp{file}.edfLocs   = (recordCmp{file}.time * epochDurationSampleRate / epochDurationSeconds);
  recordCmp{file}.tag1      = ismember(temp{file}.unionEventIndex, temp{file}.s1EventIndex);
  recordCmp{file}.tag2      = ismember(temp{file}.unionEventIndex, temp{file}.s2EventIndex);
  combined = recordCmp{file}.tag1 + recordCmp{file}.tag2;
  for ww = 1:length(combined)
      if combined > 0
          recordCmp{file}.combined(ww)  = 1;
      else
          recordCmp{file}.combined(ww) = 0;
      end
  end
  
  % exact seizure locations
  recordCmp{file}.exactSeizureStart   = zeros(length(recordCmp{file}.time), 1);
  recordCmp{file}.exactSeizureStop    = zeros(length(recordCmp{file}.time), 1);
  signal{file}                        = temp{file}.edf.signalCell{1, 2};
  
end

clear temp

% normalize the raw signal as we did in the manuscript
for i = 1:size(signal, 2)
    dataFrame(i).signal = signal{i};
    dataFrame(i).fs = 256;
    % find any blank parts of the EDF (singal == 0) and remove them from consideration for the normalization process

    % normalize signal
    [dataFrame(i).signal, dataFrame(i).rSig, dataFrame(i).modelfit, dataFrame(i).rMu] = normalizeEEG(dataFrame(i).signal, dataFrame(i).fs);

    % calculate the first derivative of the signal and normalize that
    dataFrame(i).gSignal = gradient(dataFrame(i).signal);
    [dataFrame(i).gSignal, dataFrame(i).gSig] = normalizeEEG(dataFrame(i).gSignal, dataFrame(i).fs); % save the standard deviation of the signal for later use, note this should be 1.
    
end

for i = 1:5
    recordCmp{1,i}.edfLocs = recordCmp{1,i}.edfLocs + 1024;
end

% load all of the matrices that I previously calculated and compare the output
dataFrame(1).machineOutput = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Ronde_4_exactSeizureLocations.mat');
dataFrame(1).predictorMatrix = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Ronde_4_predictorMatrix.mat');

dataFrame(2).machineOutput = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Ronde_5_exactSeizureLocations.mat');
dataFrame(2).predictorMatrix = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Ronde_5_predictorMatrix.mat');

dataFrame(3).machineOutput = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Sherlock_6_exactSeizureLocations.mat');
dataFrame(3).predictorMatrix = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Sherlock_6_predictorMatrix.mat');

dataFrame(4).machineOutput = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Tiki_4_exactSeizureLocations.mat');
dataFrame(4).predictorMatrix = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Tiki_4_predictorMatrix.mat');

dataFrame(5).machineOutput = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Vincent_4_exactSeizureLocations.mat');
dataFrame(5).predictorMatrix = load('/Volumes/cookieMonster/mj1lab_mirror/gamma2R43Q/detectSWD_output/Vincent_4_predictorMatrix.mat');

%%  calculate the event detection score using the machine

% if you run this more than one time then you must uncomment the next line to clean up the dataFrame
%dataFrame = rmfield(dataFrame, {'starts', 'ends'});

% add the machine detected events
for i = 1:length(dataFrame)
    counter = 1;
    for j = 1:length(dataFrame(i).machineOutput.output)
        %if dataFrame(i).machineOutput.output(j).classificationScore <= 0 % if this is commented out then you get all the yes and no events, if it's < 0 then you get all events considered yes by the machine.
            dataFrame(i).starts(counter) = dataFrame(i).machineOutput.output(j).SWDLOCS(1);
            dataFrame(i).ends(counter) = dataFrame(i).machineOutput.output(j).SWDLOCS(end);
            dataFrame(i).classification(counter) = dataFrame(i).machineOutput.output(j).classificationScore <= 0;
            counter = counter + 1;
        %end
    end
end

% create a list of automatedYes events and human yes to compare -- this is hard because the machine identifies events and the humans originally  scored things as epochs.. the next section of this code matches these things up.. Of note that sometimes we are dealing with epochs and sometimes we are dealing with events from here to the bottom of this code.
% if you want the results presented in the manuscript change the .combined
% stuff to .tag1 and .tag2 and then add the confusion matricies. Also, if
% the section above is set to <=0 in line 121 then you'll get the results
% for Stage 2 identified events only but if it's commented out then you are
% getting all the events.

for i = 1:size(signal, 2)
    dataFrame(i).epochsList = zeros(size(signal{i}, 1) / 256 / 4, 1);
    for zz = 1:size(dataFrame(i).epochsList, 1)
        epochEnd = (zz * 256 * 4);
        epochStart = epochEnd - 1024 + 1;
        
        % is start within epoch?
        if sum(dataFrame(i).starts >= epochStart & dataFrame(i).starts <= epochEnd) > 0
            dataFrame(i).startWithinEpoch(zz) = 1;
        else
            dataFrame(i).startWithinEpoch(zz) = 0;
        end
        
        % is end within epoch?
        if sum(dataFrame(i).ends >= epochStart & dataFrame(i).ends <= epochEnd) > 0
            dataFrame(i).endWithinEpoch(zz) = 1;
        else
            dataFrame(i).endWithinEpoch(zz) = 0;
        end  
        
        if dataFrame(i).startWithinEpoch(zz) == 1 | dataFrame(i).endWithinEpoch(zz) == 1
            dataFrame(i).automatedYes(zz) = 1;
            

        else
            dataFrame(i).automatedYes(zz) = 0;
        end
        
        % is epoch or nextdoor epochs taggeded by humans?
        if dataFrame(i).startWithinEpoch(zz) ==  1 & dataFrame(i).endWithinEpoch(zz) == 1
            if sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd) > 0
                dataFrame(i).humanYes(zz) = 1;
            else
                dataFrame(i).humanYes(zz) = 0;
            end
        elseif dataFrame(i).startWithinEpoch(zz) ==  1 & dataFrame(i).endWithinEpoch(zz) == 0
        
            if sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd) > 0 | sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart + 1024 & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd + 1024) > 0
                dataFrame(i).humanYes(zz) = 1;
            else
                dataFrame(i).humanYes(zz) = 0;
            end
        elseif dataFrame(i).startWithinEpoch(zz) ==  0 & dataFrame(i).endWithinEpoch(zz) == 1
        
            if sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd) > 0 | sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart - 1024 & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd - 1024) > 0
                dataFrame(i).humanYes(zz) = 1;
            else
                dataFrame(i).humanYes(zz) = 0;
            end
        elseif dataFrame(i).startWithinEpoch(zz) ==  0 & dataFrame(i).endWithinEpoch(zz) == 0
            if sum(recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) >= epochStart & recordCmp{i}.edfLocs(find(recordCmp{i}.combined)) <= epochEnd) > 0
                dataFrame(i).humanYes(zz) = 1;
            else
                dataFrame(i).humanYes(zz) = 0;
            end
        else
            disp('you missed one')
        end
            
        
    end
    
    %dataFrame(i).automatedYes = [dataFrame(i).automatedYes(2:end) 0];
end


% make confusion matrix

for i = 1:length(dataFrame)
    confmats(i).confmat = confusionmat(dataFrame(i).humanYes, dataFrame(i).automatedYes);
    confusionmat(dataFrame(i).humanYes, dataFrame(i).automatedYes)
end


summedconfmat = confmats(1).confmat + confmats(2).confmat + confmats(3).confmat + confmats(4).confmat + confmats(5).confmat
sum(sum(summedconfmat))


%% now calculate the predictor variables from all of the events in each category of false negatives (lower left corner) and true positives (lower right corner) for the bottom half matrix



SWDcounter = 1;
nonSWDcounter = 1;
for i = 1:length(dataFrame)
    for j = 1:size(dataFrame(i).predictorMatrix.predictorMat, 1)
        
        if dataFrame(i).classification(j) == 1
            SWDpredictors(SWDcounter, :) = dataFrame(i).predictorMatrix.predictorMat(j,:);
            SWDcounter = SWDcounter + 1;
        else
            nonSWDpredictors(nonSWDcounter, :) = dataFrame(i).predictorMatrix.predictorMat(j,:);
            nonSWDcounter = nonSWDcounter + 1;
        end
    end
end



%{
%%


% first find all the epochs that have a humanYes that don't have a
% automated yes from the combined.

%make sure you do this with the tag set to 'combined' above
figure('units', 'inch', 'pos', [3 10 20 5]);

for i = 2:size(dataFrame, 2)
    
    % these are the locations around where humans marked yes and the machine marked didn't find anything
    
    pospos = find(dataFrame(i).humanYes == 1 & dataFrame(i).automatedYes == 0);
    
    for j = 1:length(pospos)
        epochStart = (dataFrame(i).fs * 4 * pospos(j)) - 1024;
        plot(dataFrame(i).signal(epochStart - 1024: epochStart + 1024 * 2))
        hold on
        if dataFrame(i).automatedYes(pospos(j) - 1) == 1
            plot(0, 1, 'go')
        end
        if dataFrame(i).automatedYes(pospos(j)) == 1
            plot(1024, 1, 'go')
        end
        if dataFrame(i).automatedYes(pospos(j) + 1) == 1
            plot(2048, 1, 'go')
        end
                
        plot(1024, 0, 'ro')
        plot(2048, 0, 'ro')
        hold off
        [x, y] = ginput(2);
        dataFrame(i).HYMCE_start(j) = x(1) + epochStart - 512;
        dataFrame(i).HYMCE_end(j) = x(2) + epochStart - 512;
        disp(strcat(num2str(j), {' of '}, num2str(length(pospos)), {' for file '}, num2str(i)))

        
        %waitforbuttonpress()
    end
end

%%
%}

% the code commented out above is where I manually marked the events within
% the epochs that the humans found that the machine didn't.. so the matrix
% loaded from the file below has all of those timings.
load('~/Desktop/eventsMissedBytheMachine.mat')

%% now calculate the predictor variables for these events as the algorithm would have done
fs = 256;
counter = 1;
for i = 1:size(dataFrame, 2)
    for j = 1:size(dataFrame(i).HYMCE_start, 2)
        if (dataFrame(i).HYMCE_end(j) - dataFrame(i).HYMCE_start(j)) > 2
% wavelet analyses on the full signal clips
        signalclip = dataFrame(i).signal(dataFrame(i).HYMCE_start(j) -512 - 512:dataFrame(i).HYMCE_end(j)-512 + 512);
        cwtMortOut = cwt(signalclip, 'amor', fs);

        % sixHz variables
        sixHz = sum(abs(cwtMortOut(floor(fs/6.4):floor(fs/5.22),:))); %5-8hz approximately
        sixHz_mean(counter) = mean(sixHz(513:end-512));
        sixHz_std(counter) = std(sixHz(513:end-512));
        sixHz_max(counter) = max(sixHz(513:end-512));

        % inBetween variables
        inBetween = sum(abs(cwtMortOut(floor(fs/8.53):floor(fs/6.56),:))); 
        inBetween_mean(counter) = mean(inBetween(513:end-512));
        inBetween_std(counter) = std(inBetween(513:end-512));
        inBetween_max(counter) = max(inBetween(513:end-512));
        
        % calculate the harmonic  variables
        harmonic = sum(abs(cwtMortOut(floor(fs/12.8):floor(fs/8.83),:))); %15-32 hrz approx
        harmonic_mean(counter) = mean(harmonic(513:end-512));
        harmonic_std(counter) = std(harmonic(513:end-512));
        harmonic_max(counter) = max(harmonic(513:end-512));
        
        % higherFreq variables
        higherFreq = sum(abs(cwtMortOut(floor(fs/25.6):floor(fs/13.47),:))); 
        higherFreq_mean(counter) = mean(higherFreq(513:end-512));
        higherFreq_std(counter) = std(higherFreq(513:end-512));
        higherFreq_max(counter) = max(higherFreq(513:end-512));
counter = counter + 1
        end
    end

    end

    % redistribute to the featureSet object
    featureSet = table(sixHz_mean', sixHz_std', sixHz_max', inBetween_mean', inBetween_std', inBetween_max', harmonic_mean',  harmonic_std', harmonic_max', higherFreq_mean', higherFreq_std', higherFreq_max');
    featureSet.Properties.VariableNames = {'sixHz_mean', 'sixHz_std', 'sixHz_max', 'inBetween_mean', 'inBetween_std', 'inBetween_max', 'harmonic_mean', 'harmonic_std', 'harmonic_max',  'higherFreq_mean', 'higherFreq_std', 'higherFreq_max'};


    
%% now plot the data from Figure 4H

errorbar(1:12, mean(table2array(nonSWDpredictors)), std(table2array(nonSWDpredictors)), 'bo')
hold on
errorbar(1:12, mean(table2array(featureSet)), std(table2array(featureSet)), 'co')

errorbar(1:12, mean(table2array(SWDpredictors)), std(table2array(SWDpredictors)), 'ro')
plot(1:12, mean(table2array(SWDpredictors)), 'r.', 'markersize', 20)
plot(1:12, mean(table2array(featureSet)), 'c.', 'markersize', 20)
plot(1:12, mean(table2array(nonSWDpredictors)), 'b.', 'markersize', 20)

print('~/Desktop/thefinalplotofthisfuckingpaper.pdf', '-dpdf');


%% put together table for statistics

SWDpredictors.label = repmat({'SWD'}, size(SWDpredictors, 1), 1);
nonSWDpredictors.label = repmat({'nonSWD'}, size(nonSWDpredictors, 1), 1);
featureSet.label = repmat({'notFound'}, size(featureSet, 1), 1);

newTable = [SWDpredictors; nonSWDpredictors; featureSet];


%{
%% plot an example
figure
 i = 1
plot(1:length(dataFrame(i).signal),dataFrame(i).signal, 'b')
hold on
%plot(1:length(dataFrame(i).gSignal),dataFrame(i).gSignal, 'b')

plot((1:length(dataFrame(i).humanYes)) * 1024, dataFrame(i).humanYes * 3, 'om')
plot((1:length(dataFrame(i).automatedYes)) * 1024, dataFrame(i).automatedYes * 2, 'oc')
plot((dataFrame(i).starts), 4, 'og')


% from i = 1
%xlim([1.6810e+07 1.6815e+07])

% from i = 1
%xlim([2.23795e7 2.2381e7])

% from i = 1
%xlim([9.2635e+06 9.2656e+06]);

%xlim([1240064-2048 1240064+2048])


%}

