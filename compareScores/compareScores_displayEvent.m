clear all
close all
[basePath merlinPath cookieMonster] = getUserPath();

%compare Scores for Seizures or sleep on the same sets of files

% setup workspace

%basePath = '/Volumes/drWeird/Users/jesse/Dropbox/research/jonesLab';
seizurePath = '/jonesLab_data/EEG_data/RQ/seizureScoring/';
edfPath = '/jonesLab_data/EEG_data/RQ/EDFs/';
compareScoresOutputPath = '/jonesLab_data/EEG_data/RQ/compareScores/';

% prepare table to compare records

% ronde_4 = rq, ronde_5 = ganx, sherlock 4 = RR, sherlock 6 = rr and ganx, tiki_4 = rq, tiki_5 = rq + ganx, vincent_4 = RQ, vincent_10 = rq and ganx)
% the animals below are all good to compare, other animals have wierd characteristics in Rama's scoring (probably due to tag naming in Serenia) or some other weird characteristics.
files = cellstr({'Ronde_4'; 'Ronde_5'; 'Sherlock_6'; 'Tiki_4'; 'Vincent_4'});

%% import EDF signals and .csv files from comapreScores_eventID

Fs = 1024;
epochDuration = 4;
fileStart = 1; % this should be 1 most of the time

for file = fileStart:length(files)
      
  % import edf file
  temp{file}.edffilespec                                                            = char(strcat(basePath, edfPath, files(file), '.edf'));
  [temp{file}.edf.header, temp{file}.edf.signalHeader  temp{file}.edf.signalCell]   = blockEdfLoad(temp{file}.edffilespec);
  signal{file}                                                                      = temp{file}.edf.signalCell{1, 2};

  % import .csv files... Columns are epoch number, edf time of epoch num, tag1, tag2, exact start edf time, exact end edf time
  recordCmp{file} = csvread(char(strcat(basePath, compareScoresOutputPath, files(file), '.csv')), 1, 0);
  
  
end

clear temp

%% plot the signals for each animal

file = 5;

tempSignal = signal{file};
tempRecord = recordCmp{file};
eventBuffer = 100;
betweenBuffer = zeros(50, 1)';

signalMashup = zeros(length(betweenBuffer), 1)';
colorMashup = zeros(length(betweenBuffer), 1)';

    
for event = 1:length(tempRecord(:, 1))
    
    % parse out the events from each animal
    signalRange = [tempRecord(event, 5) - eventBuffer:tempRecord(event, 6) + eventBuffer];
    x = tempSignal(signalRange)';
    
    % parse the tags from each event
    tag1 = tempRecord(event, 3);
    tag2 = tempRecord(event, 4);
    
    if tag1 == 1 & tag2 == 1
        color = 0; % 0 for both tag as yes
    elseif tag1 == 1 & tag2 == 0
        color = 50; % 1 for tag1
    elseif tag1 == 0 & tag2 == 1
        color = 100; % -1 for tag2
    end

    
    % build vector for coloring by tag type
    colorVector = repmat(color, length(x), 1)';
    
    % add to mashup vector to be used for plotting
    signalMashup = [signalMashup x betweenBuffer];
    colorMashup = [colorMashup colorVector betweenBuffer];
    
end

% plot mashup for each file


% break up the mashup signal into lengths that can be plotted with the mesh function
breaks = 10000;

if length(signalMashup) < breaks
    
    fig = plot(signalMashup);
    
else
    
    % process signal for mesh plot
    signalMashupExtended = signalMashup;
    tempLength = round2(length(signalMashupExtended), breaks);
    if tempLength < length(signalMashup)
        tempLength = tempLength + breaks;
    end
    signalMashupExtended(numel(1:tempLength)) = 0;
    ncol = tempLength / breaks;
    matrixSize = [breaks, ncol];
    signalMashupMatrix = reshape(signalMashupExtended, matrixSize);
    
    % process color for mesh plot
    colorMashupExtended = colorMashup;
    tempLength = round2(length(colorMashupExtended), breaks);
    if tempLength < length(colorMashup)
        tempLength = tempLength + breaks;
    end
    colorMashupExtended(numel(1:tempLength)) = 0;
    ncol = tempLength / breaks;
    matrixSize = [breaks, ncol];
    colorMashupMatrix = reshape(colorMashupExtended, matrixSize);
    
    % process y axis for mesh plot
    tempLength = 1:breaks;
    tempLength = tempLength';
    tempY = repmat(tempLength, 1, ncol);
        
    % plot using mesh command
    colormap([1 0 0; 0 1 0; 0 0 1]) %// apply colormap, red is both, green is tag1, blue is tag2
    fig = mesh([1:size(signalMashupMatrix, 2)]', tempY, signalMashupMatrix, colorMashupMatrix, 'facecolor', 'none', 'meshstyle', 'col'); % 'edgecolor', 'k'
    %colorbar
    view(90, 85);
    
end


%% calculate FFT for 0 - 100 hz

biggestSignal  = 1:4000; % I should write some code that is dynamic to actually find out how big the biggest signal is

for file = 1:length(recordCmp)
    
    % setup signal and record objects for loop
    tempSignal = signal{file};
    tempRecord = recordCmp{file};
    eventBuffer = 1000;
    
    for event = 1:length(tempRecord(:, 1))

        % parse out the events from each animal
        signalRange = [tempRecord(event, 5) - eventBuffer:tempRecord(event, 6) + eventBuffer];
        x = tempSignal(signalRange);
        
        % add 0s to the end of the signal so they all match in length
        x2 = x;
        x2(numel(biggestSignal)) = 0;

        
        % calculate fft
        N = length(x2);
        X = fft(x2);
        X_mag = abs(X);
        
        %{
        N = length(x);
        X = fft(x);
        X_mag = abs(X);
        %}
        % establish conversions from bins to frequency axis I DONT THINK THIS IS CORRECT CURRENTLY
        f = (0:(N-1))*Fs/N;
        
        %{
        % plot the fft
        figure('units', 'inch', 'pos', [15 8 5 8]);
        plot(X_mag)
        xlim([0 100])
        title('Magnitude of DFT')
        xlabel('Frequency (Hz)')
        %}
        
        % generate index for the fft up to 100 hz. this index is dynamic based on the length of the initial event
        index = find(f <= 100);
     
        signalFFT{file}.X(:, event) = x2;
        signalFFT{file}.X_mag(:, event) = X_mag(index); 
        signalFFT{file}.freq(:, event) = f(index);
        

    end
    
    clear tempSignal tempRecord x x2 signalRange f index X_mag X N
    
end

% plot the fft events right now this only works for files with relatively few events... Perhaps change to the plotting method that we built above
file = 3;
fig = mesh([1:size(signalFFT{file}.X_mag, 2)]', signalFFT{file}.freq, signalFFT{file}.X_mag, 'edgecolor', 'k', 'facecolor', 'none', 'meshstyle', 'col');
view(90, 85);

% plot raw data
file = 1;
fig = mesh([1:size(signalFFT{file}.X, 2)]', biggestSignal, signalFFT{file}.X, 'edgecolor', 'k', 'facecolor', 'none', 'meshstyle', 'col');
view(90, 85);

%% pull a bunch of important characteristics from each of these events that we can cluster or something like that... This will give us information on the difference between events as tagged by Rama and myself as well as information on how to automatically detect events.

%% pca or other ordination of fft signals

file = 1;

pca(signalFFT{file}.X_mag)

HeatMap(signalFFT{file}.X_mag')

imagesc(zscore(signalFFT{file}.X_mag)')

%% other things

%spectrogram(tempSignal(signalRange), 'yaxis');

%imagesc(tempSignal(signalRange)');

