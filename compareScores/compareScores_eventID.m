clear all
close all
[basePath merlinPath] = getUserPath();

%compare Scores for Seizures or sleep on the same sets of files

% setup workspace

%basePath = '/Volumes/drWeird/Users/jesse/Dropbox/research/jonesLab';
seizurePath = '/jonesLab_data/EEG_data/RQ/seizureScoring/';
edfPath = '/jonesLab_data/EEG_data/RQ/EDFs/';
compareScoresOutputPath = '/jonesLab_data/EEG_data/RQ/compareScores/';

% prepare table to compare records

% ronde_4 = rq, ronde_5 = ganx, sherlock 4 = RR, sherlock 6 = rr and ganx, tiki_4 = rq, tiki_5 = rq + ganx, vincent_4 = RQ, vincent_10 = rq and ganx)
files = cellstr({'Ronde_4'; 'Ronde_5'; 'Sherlock_6'; 'Tiki_4'; 'Vincent_4'});

%%

sampleRate = 1024;
epochDuration = 4;
fileStart = 1; % this should be 1 most of the time

for file = fileStart:length(files)
  
  % filepath for each set of files
   temp{file}.edffilespec   = char(strcat(basePath, edfPath, files(file), '.edf'));
   temp{file}.s1filespec    = char(strcat(basePath, seizurePath, 'rm/', files(file), '.txt'));
   temp{file}.s2filespec 	= char(strcat(basePath, seizurePath, 'jp/', files(file), '.tag'));
  
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
  recordCmp{file}.edfLocs   = (recordCmp{file}.time * sampleRate / epochDuration);
  recordCmp{file}.tag1      = ismember(temp{file}.unionEventIndex, temp{file}.s1EventIndex);
  recordCmp{file}.tag2      = ismember(temp{file}.unionEventIndex, temp{file}.s2EventIndex);
  
  % exact seizure locations
  recordCmp{file}.exactSeizureStart   = zeros(length(recordCmp{file}.time), 1);
  recordCmp{file}.exactSeizureStop    = zeros(length(recordCmp{file}.time), 1);
  signal{file}                        = temp{file}.edf.signalCell{1, 2};
  
 
end

clear temp

%% tag exact seizure locations 

fileStart = 1;
key = 1;
loopStart = 1;
i = 1;

  
for file = fileStart:length(files)
  while i < length(recordCmp{file}.edfLocs)
      for i = loopStart:length(recordCmp{file}.edfLocs)
          tagStart = recordCmp{file}.edfLocs(i) + (sampleRate); % not sure why I have to add the sample rate here to shift the window... figure this out!
          if tagStart < 3 * sampleRate
              recordCmp{file}.exactSeizureStart(i) = 0; % place start locations in object
              recordCmp{file}.exactSeizureStop(i) = 0;
          else
              beforeWin = sampleRate * 2 - 1;
              afterWin = sampleRate * 2 - 1;
              seq = [(tagStart - beforeWin):(tagStart + afterWin)];
              fig = figure('units', 'inch', 'pos', [1 10 50 6]);
              p = plot(signal{file}(seq));
              minimum = min(signal{file}(seq)) + (min(signal{file}(seq)) * 0.1); % minimum value for the EEG seq we are looking at plus 10 percent
              maximum = max(signal{file}(seq)) + (max(signal{file}(seq)) * 0.1); % same but for max
              bottom = -1 * max(abs([minimum, maximum])); % calculates where the bottom of the rectange that highlights the active epoch should go
              top = 2 * abs(bottom); % determines the top of the rectange
              r = rectangle('Position', [sampleRate, bottom, sampleRate, top], 'EdgeColor', 'r'); % places a rectangle on the plot that indicates the orignial tag location
              set(gca, 'xlim', [0, sampleRate*3]);
              text(1100, abs(bottom) - 0.01, ['record ', num2str(i), ' of ', num2str(length(recordCmp{file}.edfLocs))]);
              if recordCmp{file}.tag1(i) == 1
                  text(1100, bottom + .01, 'tag1');
              end
              if recordCmp{file}.tag2(i) == 1
                  text(1150, bottom + .01, 'tag2');
              end
              %waitforbuttonpress;
              [a, b, key] = ginput(2); % code that allows user to input 2 clicks for start and stop
              recordCmp{file}.exactSeizureStart(i) = round(a(1,1)) + (tagStart - beforeWin); % place start locations in object
              recordCmp{file}.exactSeizureStop(i) = round(a(2,1)) + (tagStart - beforeWin); % place stop locations in object

              close all

              if key ~= 1
                  loopStart = i - 1;
                  break
              end
          end
      end
  end
end

%% export files

fileStart = 1;
for file = fileStart:length(files)
  
  tempTable = recordCmp{file};
  tempPath = char(strcat(basePath, compareScoresOutputPath, files(file), '.csv'));
  if exist(tempPath) > 0
      delete(tempPath);
      writetable(tempTable, tempPath);
  else
      writetable(tempTable, tempPath);
  end

end

%% save the file for later use

dateStamp = clock;
save([basePath, compareScoresOutputPath, 'compareScores_rm_jp_', num2str(dateStamp(2)), '_' , num2str(dateStamp(3)), '_' , num2str(dateStamp(1))], 'recordCmp');



