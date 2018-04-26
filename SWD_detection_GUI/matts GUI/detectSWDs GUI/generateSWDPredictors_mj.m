function [ features ] = generateSWDPredictors( alreadyImportedEDF, tagfilespec )
% This function generates a matrix of 19 predictors that can be used with
% the the trained classifier 'polynomialSWDSVMClassifier.m'.
%
%   This function takes either 1 input, an already imported EDF (using
%   read_EDF_mj), or 2. In the instance with 1 inputs, it generates a
%   matrix of predictor variables. If a 2rd input is supplied (tagfilespec, a filepath for a
%   svarog or sirenia sleep pro file tagged for seizure locations) then it
%   outputs a matrix of predictors with SWD labels and one other thing in the first
%   column two columns. This script was designed for use in developing and SVM for
%   predicting SWD epochs in EEG and predicting SWDs in novel data.
%
%      
%
%   Example from Ronde_4.edf, i = 666 is a good SWD epoch.
%
%   [basePath merlinPath cookieMonster] = getUserPath();
%
%   % load EDF
%   edffilespec = strcat(basePath, '/jonesLab_tempData/sleep_and_seizures/Ronde_4.edf');
%   alreadyImportedEDF = read_EDF_mj(edffilespec);
%
%   % optionally load tagfilespec
%   tagfilespec = strcat(basePath, '/jonesLab_tempData/sleep_and_seizures/Ronde_4.tag');
%
%   % generate predictor matrix
%   SWDPredMatrix = generateSWDPredictors( EDF, tagfilespec );
%
%   JP 2016   


% variable argument defaults
if nargin > 2
    error('Function takes 1 or 2 arguments.')
end
switch nargin
    case 1
        includeSWDtags = 0;
    case 2
        includeSWDtags = 1;
end


% set paths
[basePath, ~, ~] = getUserPath();

% Abbreviated analysis parameters based on the current file
fs = alreadyImportedEDF.D.edf.fs;                           % sample frequency
epochdur = alreadyImportedEDF.D.edf.epochdur;
nepochs = alreadyImportedEDF.D.edf.nsecs./alreadyImportedEDF.D.edf.epochdur;

% label the signal from the edf
signal = normalizeEEG(alreadyImportedEDF.D.edf.eegData);    % Raw signal
alreadyImportedEDF.D.edf.eegData = signal;                  % add back to raw signal for later use in power stuff.

% tags for signal definitions
FR = 1;     % F for front, left and right for the second letter. check order or left and right
%FL = 2;     % this seems to be the strongest SWD signal for our animals
%PR = 3;    % parietal right -- not sure if I got the order of right and left correct or not.
%PL = 4;
EMG = 5;
pixel = 6;

% SWD Template
SWDtemplatefilespec = strcat(basePath, '/jonesLab_git/jonesLab_code/Projects/sleep_and_seizures/EEG_Analysis/seizureIdentification/characterizeSWDs/SWDTemplate.csv');
SWDTemplate = csvread(SWDtemplatefilespec); % this template was generated from 'designSWDtemplate.m'

%{
sleepLabel = alreadyImportedSleepScores.score;
NREM = strcmp(sleepLabel, 'Non REM');
REM = strcmp(sleepLabel, 'REM');
Wake = strcmp(sleepLabel, 'Wake');
%}
% import SWD tags if an agrument was provided
if includeSWDtags == 1
    % import seizure score
    S = read_svarog_tagfile(tagfilespec); % be sure that svarog tagfiles have start and stop tag at beginning and end.

    lasttaggedepoch = S.tagdata.epochtime(end)./epochdur;
    S.tagdata.isSWD(lasttaggedepoch+1:nepochs) = 0;
    S.tagdata.epochtime(lasttaggedepoch+1:nepochs) = S.tagdata.epochtime(end)+epochdur.*(nepochs-lasttaggedepoch);

    % sleep and SWD indicies
    SWDLabs = logical(S.tagdata.isSWD);
    possibleSWDLabsInd = S.tagdata.tagepoch(strcmp(S.tagdata.tagname, 'Possible SWD'));
    possibleSWDLabs = zeros(length(SWDLabs), 1);
    for i = 1:length(possibleSWDLabsInd)
        j = possibleSWDLabsInd(i);
        possibleSWDLabs(j) = 1;
    end

end


% Generate predictor variables for SVM

%% ---- wavelet analysis
thresh = 5;

%{
redo = 1;
while redo == 1
%}
wb = waitbar(0,'1/5 - Computing wavelet decompositions...');
    for i = 1:nepochs 
            waitbar(i/nepochs,wb)

        % setups for the signal
        nepochpts = fs * epochdur;
        epochEnd = i * nepochpts;
        epochStart = epochEnd - nepochpts + 1;
        signalI = signal(epochStart:epochEnd, :);


        dt = 1/256; % sample rate = 256
        %t = linspace(0,epochdur,nepochpts);
        cwt = cwtft({signalI(:,FR)', dt});
        %cwt = cwtft({signalI(:,FR)', dt}, 'plot'); % this shows good signal for i = 666 SWD in Ronde_4

        %contour(t, cwt.frequencies, real(cwt.cfs), '.');

        % the contour plot for the SWD in i = 666 for Ronde_4 shows a peak in the block between 5 and 7, which is represented by row 10 @ ~5.8 hz
        SWDWaveletSignal = abs(cwt.cfs(10,:));
        
        
        sumSWDWaveletSignal(i) = sum(SWDWaveletSignal);
        sumThreshCrossingsSWDWaveletSignal(i) = sum(SWDWaveletSignal > thresh);
        threshCrossingsSWDWaveletSignal(i) = sum(diff(SWDWaveletSignal > thresh) > 0);

        %{
        % use findpeaks to find the number of threshold crossings
        warning('off', 'signal:findpeaks:largeMinPeakHeight');
        [pks, locs] = findpeaks(SWDWaveletSignal, 'MinPeakHeight', MPH);
        nWaveletPKS(i) = length(pks);

        if length(locs) == 0
            meanWaveletPKS(i) = 0;
            meanWaveletLOCS(i) = 0;
            varWaveletLOCS(i) = 0;
        elseif length(locs) == 1
            meanWaveletPKS(i) = mean(pks);
            meanWaveletLOCS(i) = 0;
            varWaveletLOCS(i) = 0;  
        else
            meanWaveletPKS(i) = mean(pks);
            meanWaveletLOCS(i) = mean(diff(locs));
            varWaveletLOCS(i) = var(diff(locs));
        end
        %}
    end
close(wb)
%{
    check = confusionmat(SWDLabs, threshCrossingsSWDWaveletSignal > 0);
    threshWavlet = thresh;
    thresh = thresh * 0.95; % perhaps I should set this threshold differently instead of manually
    redo = check(2, 1) / check(2, 2) > 0; % keep going unless you find all but 5 percent of the SWDs
end
%}

%% ---- use the differential of the signal to identify SWDs, this seems to be a pretty good predictor
MPH = .8574; % threshold set after starting at 1 and using the while loop below.

%{
redo = 1;
while redo == 1
%}
wb = waitbar(0,'2/5 - Differentiating signals...');
    for i = 1:nepochs 
        waitbar(i/nepochs,wb)
        % setups for the signal
        nepochpts = fs * epochdur;
        epochEnd = i * nepochpts;
        epochStart = epochEnd - nepochpts + 1;
        signalI = signal(epochStart:epochEnd, :);

        diffSignalI = diff(signalI(:,FR));
        warning('off', 'signal:findpeaks:largeMinPeakHeight');
        [pks, locs] = findpeaks(diffSignalI, 'MinPeakHeight', MPH);
        nDiffThreshPKS(i) = length(pks);

        if length(locs) == 0
            meanDiffThreshPKS(i) = 0;
            meanDiffThreshLOCS(i) = 0;
            varDiffThreshLOCS(i) = 0;
        elseif length(locs) == 1
            meanDiffThreshPKS(i) = mean(pks);
            meanDiffThreshLOCS(i) = 0;
            varDiffThreshLOCS(i) = 0;  
        else
            meanDiffThreshPKS(i) = mean(pks);
            meanDiffThreshLOCS(i) = mean(diff(locs));
            varDiffThreshLOCS(i) = var(diff(locs));
        end
    end
close(wb)
%{
    check = confusionmat(SWDLabs, logical(nDiffThreshPKS));
    MPHdiff = MPH;
    MPH = MPH * 0.95; % perhaps I should set this threshold differently instead of manually
    redo = check(2, 1) / check(2, 2) > 0; % keep going unless you find all but 5 percent of the SWDs
end
%}
    
%% ----- template matching for SWD
MPH = 1.0267; % set initial threshold at 2 and allowed to go down to 1.0267 with redo code below

%{
redo = 1;
while redo == 1
%}
wb = waitbar(0,'3/5 - Template Matching...');
    for i = 1:nepochs
        waitbar(i/nepochs,wb)

        % setups for the signal
        nepochpts = fs * epochdur;
        epochEnd = i * nepochpts;
        epochStart = epochEnd - nepochpts + 1;
        signalI = signal(epochStart:epochEnd, :);

        meanSignalI = mean(signalI(:, 1:4), 2);

        corr = xcorr(meanSignalI, SWDTemplate);
        warning('off', 'signal:findpeaks:largeMinPeakHeight');
        [PKS, LOCS] = findpeaks(corr, 'MinPeakHeight', MPH);

        % number of peaks
        nTemplatePKS(i) = length(PKS);


        % mean and variance of the distance between peaks
        if length(LOCS) == 0
            meanTemplatePKS(i) = 0;
            meanTemplateLOCS(i) = 0;
            varTemplateLOCS(i) = 0;
        elseif length(LOCS) == 1
            meanTemplatePKS(i) = mean(PKS);
            meanTemplateLOCS(i) = 0;
            varTemplateLOCS(i) = 0;
        else
            meanTemplatePKS(i) = mean(PKS);
            meanTemplateLOCS(i) = mean(diff(LOCS));
            varTemplateLOCS(i) = var(diff(LOCS));
        end
    end
close(wb)
%{
    check = confusionmat(SWDLabs, logical(nTemplatePKS));
    MPHtemplate = MPH;
    MPH = MPH * 0.95; % perhaps I should set this threshold differently instead of manually
    redo = check(2, 1) / check(2, 2) > 0; % keep going unless you find all but 5 percent of the SWDs
end
%}
    
    
%{
    %% coherence between signals -- should SWDs have high coherence??

for i = 1:nepochs
    
    % setups for the signal
    nepochpts = fs * epochdur;
    epochEnd = i * nepochpts;
    epochStart = epochEnd - nepochpts + 1;
    signalI = signal(epochStart:epochEnd, :);
   
    % correlations
    %correlations12(i) = corr(signalI(:, 1), signalI(:, 2));
        correlations12(i) = mscohere(signalI(:, 1), signalI(:, 2));

    correlations23(i) = corr(signalI(:, 2), signalI(:, 3));
    correlations34(i) = corr(signalI(:, 3), signalI(:, 4));

    
end

plot3(correlations12, correlations23, correlations34, '.')
hold on;
plot3(correlations12(SWDLabs), correlations23(SWDLabs), correlations34(SWDLabs), 'ro')
hold off;
    
%}

%% more predictors that don't have a threshold
wb = waitbar(0,'4/5 - Additional Predictors...');
for i = 1:nepochs 
	waitbar(i/nepochs,wb)
    % setups for the signal
    nepochpts = fs * epochdur;
    epochEnd = i * nepochpts;
    epochStart = epochEnd - nepochpts + 1;
    signalI = signal(epochStart:epochEnd, :);


    %% ----- variance estimates of the signal mean, this seems to be a decent predictor as well.
    meanSignalI = mean(signalI(:, 1:4), 2);
    var_meanSignalI(i) = var(meanSignalI);
    mean_meanSignalI(i) = mean(abs(meanSignalI));

    %{
    subplot(2, 1, 1)
        plot(theMatrix.varOfMeanSignal, theMatrix.SWDs, '.') %
    subplot(2, 1, 2)
        histogram(theMatrix.varOfMeanSignal)
    %}


    %% ----- Sum of Pixels that signify movement
    sumPixel(i) = sum(signalI(:, pixel));

    %{
    subplot(2, 1, 1)
        plot(theMatrix.sumPixels, theMatrix.SWDs, '.') % possibly helpful
    subplot(2, 1, 2)
        histogram(theMatrix.sumPixels)
    %}

    %% ----- Sum of EMG signal
    sumEMG(i) = sum(abs(signalI(:, EMG)));

    %{
    subplot(2, 1, 1)
        plot(theMatrix.sumEMG, theMatrix.SWDs, '.') % potentially helpful
    subplot(2, 1, 2)
        histogram(theMatrix.sumEMG)
    %}

end
close(wb)

    %% ----- get all edf power stuff (some calculations are repeated here from earlier in the script)
    powerStuff = get_all_edf_power_jp(alreadyImportedEDF);
    bandPowers = powerStuff.D.edf.prepostprocdata.dump.EDFstats;
    clear powerStuff

    avg_Mag = cell2mat(squeeze(bandPowers.avg_Mag(1,:,:))');
    max_Mag = cell2mat(squeeze(bandPowers.max_Mag(1,:,:))');
    %{
    delta   = [0.5 4];                  % Freq ranges for bands to analyze
    SWD     = [5.5 9.5];
    SWDHarm = [11 18];                
    theta   = [6 9];   
    sigma   = [10 14];   
    gamma   = [25 100];
    %}
    
    %
    avgDelta = avg_Mag(:,1);
    avgTheta = avg_Mag(:,4);
    avgSigma = avg_Mag(:,5);
    avgGamma = avg_Mag(:,6);
    
    deltaGammaRatio_avgMag = avgDelta ./ avgGamma;
    sigmaGammaProduct_avgMag = avgSigma .* avgGamma;
    
    maxDelta = max_Mag(:,1);
    maxTheta = max_Mag(:,4);
    maxSigma = max_Mag(:,5);
    maxGamma = max_Mag(:,6);
    
    deltaGammaRatio_maxMag = avgDelta ./ avgGamma;
    sigmaGammaProduct_maxMag = avgSigma .* avgGamma;
    
    %% -----  concatonate matrix
    features = [sumSWDWaveletSignal; sumThreshCrossingsSWDWaveletSignal; threshCrossingsSWDWaveletSignal; nDiffThreshPKS; meanDiffThreshPKS; meanDiffThreshLOCS; varDiffThreshLOCS; var_meanSignalI; mean_meanSignalI; sumPixel; sumEMG; nTemplatePKS; meanTemplatePKS; meanTemplateLOCS; varTemplateLOCS]';
    
    if includeSWDtags == 1
        features = horzcat(SWDLabs(1:end-1)', possibleSWDLabs(1:end-1), features(1:end-1,:), avgDelta, avgTheta, avgSigma, avgGamma, deltaGammaRatio_avgMag, sigmaGammaProduct_avgMag, maxDelta, maxTheta, maxSigma, maxGamma, deltaGammaRatio_maxMag, sigmaGammaProduct_maxMag);
        features = array2table(features, 'VariableNames', {'SWDs', 'PossibleSWDLabs',  'sumSWDWaveletSignal', 'sumThreshCrossingsSWDWaveletSignal', 'threshCrossingsSWDWaveletSignal', 'nDiffThreshPKS', 'meanDiffThreshPKS', 'meanDiffThreshLOCS', 'varDiffThreshLOCS', 'varOfMeanSignal', 'meanOfMeanSignal', 'sumPixels', 'sumEMG', 'nSWDTemplatePKS', 'meanSWDTemplatePKS', 'meanSWDTemplateLOCS', 'varSWDTemplateLOCS', 'avgDelta', 'avgTheta', 'avgSigma', 'avgGamma', 'deltaGammaRatio_avgMag', 'sigmaGammaProduct_avgMag', 'maxDelta', 'maxTheta', 'maxSigma', 'maxGamma', 'deltaGammaRatio_maxMag', 'sigmaGammaProduct_maxMag'});
    elseif includeSWDtags == 0
        features = horzcat(features(1:end-1,:), avgDelta, avgTheta, avgSigma, avgGamma, deltaGammaRatio_avgMag, sigmaGammaProduct_avgMag, maxDelta, maxTheta, maxSigma, maxGamma, deltaGammaRatio_maxMag, sigmaGammaProduct_maxMag);
        features = array2table(features, 'VariableNames', {'sumSWDWaveletSignal', 'sumThreshCrossingsSWDWaveletSignal', 'threshCrossingsSWDWaveletSignal', 'nDiffThreshPKS', 'meanDiffThreshPKS', 'meanDiffThreshLOCS', 'varDiffThreshLOCS', 'varOfMeanSignal', 'meanOfMeanSignal', 'sumPixels', 'sumEMG', 'nSWDTemplatePKS', 'meanSWDTemplatePKS', 'meanSWDTemplateLOCS', 'varSWDTemplateLOCS', 'avgDelta', 'avgTheta', 'avgSigma', 'avgGamma', 'deltaGammaRatio_avgMag', 'sigmaGammaProduct_avgMag', 'maxDelta', 'maxTheta', 'maxSigma', 'maxGamma', 'deltaGammaRatio_maxMag', 'sigmaGammaProduct_maxMag'});
    end
    
end

%% extra code and plotting

% extra plotting for power stuff
    %{
    delta   = [0.5 4];                  % Freq ranges for bands to analyze
    SWD     = [5.5 9.5];
    SWDHarm = [11 18];                
    theta   = [6 9];   
    sigma   = [10 14];   
    gamma   = [25 100];
    %}
    
    %{
    % average magnitude from power analysis measures
    subplot(2, 1, 1)
        plot(theMatrix.avg_Mag_delta, theMatrix.SWDs, '.') % helpful potentially
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_delta)

    subplot(2, 1, 1)
        plot(theMatrix.avg_Mag_SWD, theMatrix.SWDs, '.') % this is surprisingly unhelpful metric.. why?
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_SWD)

    subplot(2, 1, 1)
        plot(theMatrix.avg_Mag_SWDHarm, theMatrix.SWDs, '.') % this is surprisingly unhelpful metric.. why?
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_SWDHarm)

    subplot(2, 1, 1) % (13)
        plot(theMatrix.avg_Mag_theta, theMatrix.SWDs, '.') % unsure
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_theta)

    subplot(2, 1, 1)
        plot(theMatrix.avg_Mag_sigma, theMatrix.SWDs, '.') % this one might be helpful
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_sigma)

    subplot(2, 1, 1)
        plot(theMatrix.avg_Mag_gamma, theMatrix.SWDs, '.')
    subplot(2, 1, 2)
        histogram(theMatrix.avg_Mag_gamma)

    % max_Mag from power analysis
    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_delta, theMatrix.SWDs, '.') % maybe helpful
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_delta)

    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_SWD, theMatrix.SWDs, '.') % this is surprisingly unhelpful metric.. why?
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_SWD)

    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_SWDHarm, theMatrix.SWDs, '.') % this is surprisingly unhelpful metric.. why?
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_SWDHarm)

    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_theta, theMatrix.SWDs, '.') % probably not helpful
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_theta)

    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_sigma, theMatrix.SWDs, '.') % probably not helpful
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_sigma)

    subplot(2, 1, 1)
        plot(theMatrix.max_Mag_gamma, theMatrix.SWDs, '.') % maybe helpful
    subplot(2, 1, 2)
        histogram(theMatrix.max_Mag_gamma)
    %%
    plot3(theMatrix.avg_Mag_delta, theMatrix.avg_Mag_gamma, theMatrix.SWDs, '.')
    set(gca, 'xsc', 'log', 'ysc', 'log');
    %% ----- ****
    plot3(theMatrix.avg_Mag_delta ./ theMatrix.avg_Mag_gamma, theMatrix.avg_Mag_sigma .* theMatrix.avg_Mag_gamma, theMatrix.SWDs, '.')
    set(gca, 'xsc', 'log', 'ysc', 'log');
    xlabel('delta/gamma');

    %}

% frequency decomposition metrics (Modify to scan FFT for harmonics??)
%{
bin_vals = 0 : nepochpts-1;              	% Freq bins
fax_Hz = bin_vals*fs/nepochpts;             % Freq Hz
nepochpts_2 = ceil(nepochpts/2);    
X_mags = abs(fft(signalI(:, FR)));          % Simple Power Spectrum
X_mags = X_mags(1:nepochpts_2);
f = fax_Hz(1:nepochpts_2);
[~, ind] = max(X_mags); % does this need to be standardized to be able to compare epochs?
freq_at_max(i) = f(ind);

%}


% confusion matrix to check for various performance metrics
%confusionmat(logical(SWDPredMatrix.SWDs), logical(SWDPredMatrix.nWaveletPKS))



    