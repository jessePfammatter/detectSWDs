%Load this shit

% ----- import the namefile and load the aggregated data
namefilespec = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/sleep_and_seizures_RQ_namefile.txt');
F = ParseNameFile(namefilespec); % makes the object 'F'
load(strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/detectSWD_output/aggregated_RQ_Data.mat'));

% ----- 

colors = colormap(redgreenblue());
%% Figure 5

for i = 1%1:length(F.itemparams)
    for l = 3%1:length(F.itemparams{1, i}.recordDay)
        j = str2num(F.itemparams{1, i}.recordDay{1, l});
        
        % ----- PREPPING DATA FOR PLOTS
        
        % set up some variables from the namefile
        animalID = F.itemparams{1, i}.animalID{1};
        edffilespec = animal(i).day(j).edffilespec;
        [path, name, ext] = fileparts(edffilespec);
        
        % edf for spike triggered averages
        alreadyImportedEDF = read_EDF_mj(edffilespec);
        edf = alreadyImportedEDF.D.edf.signalMat(:, 2);
        fs = alreadyImportedEDF.D.edf.fs;
        
        % time variables for use in plotting
        hourInSeconds = 1:3600;
        totalHours = 1:24;
        twentyfourhours = 86400;
        plotPositions = totalHours - .5;
        
        % labely stuff for use in plotting
        hoursLabels = {'8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '1', '2', '3', '4', '5', '6', '7'};
        hoursTicks = hourInSeconds(end) * (1:24);
        
        % does sleep score and a analyzed file for swds?
        if ~strcmp(F.itemparams{1, i}.sleepScore{1, l}, 'NA') && str2num(F.itemparams{1, i}.SWDdetectOK{1, l}) == 1
        
            % set up sleep stage stuff
            sleepScoresText = animal(i).day(j).sleepScore.score;
            wakeInd = find(strcmp(sleepScoresText, 'Wake'));
            REMInd = find(strcmp(sleepScoresText, 'REM'));
            nREMInd = find(strcmp(sleepScoresText, 'Non REM'));
            artInd = find(strcmp(sleepScoresText, 'Artifact'));
            sleepScores = zeros(length(sleepScoresText), 1);
            sleepScores(wakeInd) = 1;
            sleepScores(REMInd) = 3;
            sleepScores(nREMInd) = 2;
            sleepScores(artInd) = NaN;
            time = animal(i).day(j).sleepScore.time;         
            
            % setup up some stuff for plotting SWDs
            SWDs = animal(i).day(j).SWDs;

            for ii = 1:size(SWDs, 2)
                SWDStartLOCS(ii) = SWDs(ii).SWDLOCS(1); 
                
            end
           
            %{
            % filter out all the ones based on the humanScored index that were 100 percent not SWDs as manually decided.
            SWDStartLOCS = SWDStartLOCS(find(animal(i).day(j).classificationResponse == 1));
            nonSWDStartLOCS = SWDStartLOCS(find(animal(i).day(j).classificationResponse == 2));
            %}
            
            % now expand this record to a mark for every second in which an SWD was triggered.
            SWDStartLOCS = floor(SWDStartLOCS/256); % convert to seconds instead of sample points
            everySecond = 1:time(end); % make a vector that is 24 hours worth in seconds
            SWDsecond = zeros(1, size(everySecond, 2))'; % initialize ones for SWD presence at these seconds
            for ii = 1:length(SWDStartLOCS)
                SWDsecond(SWDStartLOCS(ii)) = 1; % now add a SWD where they occur in seconds
            end
            
            % drop pairs of SWDs that start in the same second..
            if sum(SWDsecond) ~= length(SWDStartLOCS)
               diffSWDStartLOCS = [ 10, diff(SWDStartLOCS)]; % 10 is just a placeholder bigger than the next threshold
               keepInd = find(diffSWDStartLOCS);
               SWDStartLOCS = SWDStartLOCS(keepInd);
            else
                keepInd = 1:sum(SWDsecond);
            end
               
            
            % expand the sleep scores from epochs to an every second record. everySecondWake, Sleep and nREM
            sleepScoresExpanded = [];
            for ii = 1:length(sleepScores) - 1
                sleepScoresExpanded((ii * 4)) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 1) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 2) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 3) = sleepScores(ii);
            end
               
            % find the intersect of SWDs occuring at each sleep score and produce an index
            wiswd = sleepScoresExpanded == 1 & SWDsecond' == 1;
            riswd = sleepScoresExpanded == 3 & SWDsecond' == 1;
            nriswd = sleepScoresExpanded == 2 & SWDsecond' == 1;
            
            % make a matrix of the SWD characteristics to plot intensity -- perhaps amplitude and duration?
            SWDInd = find(SWDsecond);
            critsArray = zeros(1, size(everySecond, 2))'; % this is a zeros index so that 
            classificationScore = animal(i).day(j).classificationScore;
            classificationScore = classificationScore(keepInd);
            classificationResponse = animal(i).day(j).classificationResponse;
            classificationResponse = classificationResponse(keepInd);
            
            
            % add negative value for events that are NOT SWDs
            for ii = 1:length(keepInd)
                if classificationResponse(ii) == 1
                    classificationScore(ii) = -classificationScore(ii);
                end
            end
            
            critsArray(SWDInd) = classificationScore;

           
                
            
            % SWD AVERAGES PER HOUR
            yesYes = SWDsecond & (critsArray > 0); % you can change this threshold to a more negative number to include more events and a more positive number to include fewer events.
            
            WAKE_SWDsecond = yesYes' == 1 & wiswd == 1;
            REM_SWDsecond = yesYes' == 1 & riswd == 1;
            nREM_SWDsecond = yesYes' == 1 & nriswd == 1;
            for ii = totalHours
                hourInd = hourInSeconds + (hourInSeconds(end) * ii) - hourInSeconds(end);
                seizuresPerHour(ii) = sum(yesYes(hourInd));
                WAKE_seizuresPerHour(ii) = sum(WAKE_SWDsecond(hourInd));
                REM_seizuresPerHour(ii) = sum(REM_SWDsecond(hourInd));
                nREM_seizuresPerHour(ii) = sum(nREM_SWDsecond(hourInd));
            end
            
            % SLEEP STAGE AVERAGES PER HOUR
            for ii = totalHours
                hourInd = hourInSeconds + (hourInSeconds(end) * ii) - hourInSeconds(end);
                WAKE_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 1);
                percentage_WAKE_perHour(ii) = WAKE_perHour(ii);
                REM_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 3);
                percentage_REM_perHour(ii) = REM_perHour(ii);
                nREM_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 2);
                percentage_nREM_perHour(ii) = nREM_perHour(ii); 
            end
            % convert to percentages
            percentage_WAKE_perHour = percentage_WAKE_perHour / hourInSeconds(end);
            percentage_REM_perHour = percentage_REM_perHour / hourInSeconds(end);
            percentage_nREM_perHour = percentage_nREM_perHour / hourInSeconds(end);
            

            % ----- PLOTTING
            
            % open a figure
            figure('units', 'inch', 'pos', [10 10 15 15]);
            % add the main title
            h = suptitle({name});
            set(h, 'Interpreter', 'none');
            
            % PLOT 1, create a plot window for sleep stages over time
            subplot(9, 4, 1:8)
                % add a patch for lights on/lights off
                patch([hoursTicks(11), hoursTicks(23), hoursTicks(23), hoursTicks(11)], [1, 1, 3, 3], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                hold on
                text(hoursTicks(5), 3.1, 'LIGHTS ON');

                % plot the sleep scores
                plot(time, sleepScores, 'k')
                
                % add circles for sleep stages
                plot(time(wakeInd), sleepScores(wakeInd), '.c', 'MarkerSize', 10)
                plot(time(REMInd), sleepScores(REMInd), '.m', 'MarkerSize', 10)
                plot(time(nREMInd), sleepScores(nREMInd), '.g', 'MarkerSize', 10)
                
                % limits and labels
                ylim([1, 3])
                ylabel('1 = Wake, 2 = nREM, 3 = REM', 'interpreter', 'none');
                xlim([0 twentyfourhours])
                xlabel('Time of Day (24-hr)')
                xticks(hoursTicks);
                xticklabels(hoursLabels);
                
                
                % add marks for when drugs are added and what the drug was
                if ~strcmp(F.itemparams{1, i}.drugTrt{1, l}, 'none')

                    injTime1 = F.itemparams{1, i}.injTimeEDF1{1, l};
                    injTime1Num = str2num(injTime1);
                    injTime2 = F.itemparams{1, i}.injTimeEDF2{1, l};
                    injTime2Num = str2num(injTime2);

                    % add text for what drug was given
                    if ~isempty(injTime1Num)
                        %plot(injTime1Num, 1, 'b.', 'MarkerSize', 40);
                        line(injTime1Num * [1, 1], get(gca, 'ylim'), 'LineWidth', 3);
                        text(injTime1Num, 3.1, strcat(F.itemparams{1, i}.drugTrt{1, l}, {' '}, F.itemparams{1, i}.drugConc{1, l}), 'Color', 'b');
                    else
                        text(55000, 3.1, strcat(F.itemparams{1, i}.drugTrt{1, l}, {' '}, F.itemparams{1, i}.drugConc{1, l}, {' UNKNOWN DELIVERY TIMES!'}));
                    end

                    if ~isempty(injTime2Num)
                        line(injTime2Num * [1, 1], get(gca, 'ylim'), 'LineWidth', 3);
                    end 
                end
               
                hold off
               
                
            % PLOT 2, create a plot window for SWDs triggered over time
            subplot(9, 4, 9:16)
                
                % add a patch for lights on/lights off
                classifierRange = 2; % this is the range of the classifier score
                patch([hoursTicks(11), hoursTicks(23), hoursTicks(23), hoursTicks(11)], [-classifierRange, -classifierRange, classifierRange, classifierRange], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                hold on

                
                % plot the SWDs with some criteria for intesnsity (right now critsArray)
                plot(everySecond, critsArray, 'r', 'LineWidth', 0.1) % SWDsecon for a verticle line at every SWD
                

                % marks for the sleep stages sleep stages -- probably could get rid of these since I made plot 4 below
                plot(everySecond(wiswd), critsArray(wiswd), '.c', 'MarkerSize', 10)
                plot(everySecond(riswd), critsArray(riswd), '.m', 'MarkerSize', 10)
                plot(everySecond(nriswd), critsArray(nriswd), '.g' , 'MarkerSize', 10)
                
                % labels and limits
                xlabel('Time of Day (24-hr)')
                ylabel('Dist. to Nearest SV (>0 = SWD)')
                ylim([-classifierRange, classifierRange]);
                xlim([0 twentyfourhours])
                xticks(hoursTicks);
                xticklabels(hoursLabels);
               
                
                % add marks for when drugs are added and what the drug was
                if ~strcmp(F.itemparams{1, i}.drugTrt{1, l}, 'none')

                    injTime1 = F.itemparams{1, i}.injTimeEDF1{1, l};
                    injTime1Num = str2num(injTime1);
                    injTime2 = F.itemparams{1, i}.injTimeEDF2{1, l};
                    injTime2Num = str2num(injTime2);

                    % add text for what drug was given
                    if ~isempty(injTime1Num)
                        %plot(injTime1Num, 1, 'b.', 'MarkerSize', 40);
                        line(injTime1Num * [1, 1], get(gca, 'ylim'), 'LineWidth', 3);
                        text(injTime1Num, 26.5, strcat(F.itemparams{1, i}.drugTrt{1, l}, {' '}, F.itemparams{1, i}.drugConc{1, l}), 'Color', 'b');
                    else
                        text(55000, 26.5, strcat(F.itemparams{1, i}.drugTrt{1, l}, {' '}, F.itemparams{1, i}.drugConc{1, l}, {' UNKNOWN DELIVERY TIMES!'}));
                    end

                    if ~isempty(injTime2Num)
                        line(injTime2Num * [1, 1], get(gca, 'ylim'), 'LineWidth', 3);
                    end 
                end
                
            % PLOT 3, plot the moving average of SWDs   
            subplot(9, 4, 17:20)
                
                % add a patch for lights on/lights off
                patch([totalHours(11), totalHours(23), totalHours(23), totalHours(11)], [0, 0, 200, 200], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                hold on
                
                % plot the SWDs
                plot(plotPositions, seizuresPerHour, 'r.-', 'MarkerSize', 40);                
                plot(plotPositions, WAKE_seizuresPerHour, 'c.-', 'MarkerSize', 20);
                plot(plotPositions, REM_seizuresPerHour, 'm.-', 'MarkerSize', 20);
                plot(plotPositions, nREM_seizuresPerHour, 'g.-', 'MarkerSize', 20);
                
                % add labels and limits
                if max(seizuresPerHour) > 15
                    ylimSWDs = max(seizuresPerHour);
                else
                    ylimSWDs = 15;
                end
                ylim([0, ylimSWDs])
                ylabel('SWDs/Hour', 'interpreter', 'none')
                xlim([0, 24]);
                xlabel('Time of Day (24-hr)');
                xticks(totalHours);
                xticklabels(hoursLabels);
                hold off
            
 
                
                % spike triggered average stuff
                
                edf = edf(1:twentyfourhours * 256);

               
                
                critsArray_above0 = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) > 0 
                        critsArray_above0(ii) = 1;
                    else
                        critsArray_above0(ii) = 0;
                    end 
                end
                newInd = find(critsArray_above0) * 256;
                spikeArray_above0 = zeros(length(edf), 1);
                spikeArray_above0(newInd) = 1;
                 
                critsArray_0to1 = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) > 0 && critsArray(ii) <= .5
                        critsArray_0to1(ii) = 1;
                    else
                        critsArray_0to1(ii) = 0;
                    end 
                end
                newInd = find(critsArray_0to1) * 256;
                spikeArray_0to1 = zeros(length(edf), 1);
                spikeArray_0to1(newInd) = 1;
                
                critsArray_1to5 = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) > .5 && critsArray(ii) <= 1
                        critsArray_1to5(ii) = 1;
                    else
                        critsArray_1to5(ii) = 0;
                    end 
                end     
                newInd = find(critsArray_1to5) * 256;
                spikeArray_1to5 = zeros(length(edf), 1);
                spikeArray_1to5(newInd) = 1;
                
                critsArray_5toMax = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) > 1
                        critsArray_5toMax(ii) = 1;
                    else
                        critsArray_5toMax(ii) = 0;
                    end 
                end      
                newInd = find(critsArray_5toMax) * 256;
                spikeArray_5toMax = zeros(length(edf), 1);
                spikeArray_5toMax(newInd) = 1;
                
                critsArray_0toneg1 = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) < 0 && critsArray(ii) >= -.5
                        critsArray_0toneg1(ii) = 1;
                    else
                        critsArray_0toneg1(ii) = 0;
                    end 
                end       
                newInd = find(critsArray_0toneg1) * 256;
                spikeArray_0toneg1 = zeros(length(edf), 1);
                spikeArray_0toneg1(newInd) = 1;

                critsArray_neg1toneg5 = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) <= -.5 && critsArray(ii) > -1
                        critsArray_neg1toneg5(ii) = 1;
                    else
                        critsArray_neg1toneg5(ii) = 0;
                    end 
                end       
                newInd = find(critsArray_neg1toneg5) * 256;
                spikeArray_neg1toneg5 = zeros(length(edf), 1);
                spikeArray_neg1toneg5(newInd) = 1;
                
                critsArray_neg5toMin = critsArray;
                for ii = 1:length(critsArray)
                    if critsArray(ii) <= -1
                        critsArray_neg5toMin(ii) = 1;
                    else
                        critsArray_neg5toMin(ii) = 0;
                    end 
                end       
                newInd = find(critsArray_neg5toMin) * 256;
                spikeArray_neg5toMin = zeros(length(edf), 1);
                spikeArray_neg5toMin(newInd) = 1;
                
                % STA for sleep stuff
                
                % find edfpower
                powerOutput = singleChannelPower(edf, 256, 4);
                deltaPower = powerOutput.avg_Mag(:,1);
                deltaPower = deltaPower/max(deltaPower);
                thetaPower = powerOutput.avg_Mag(:,4);
                thetaPower = thetaPower/max(thetaPower);
                sigmaPower = powerOutput.avg_Mag(:,5);
                sigmaPower = sigmaPower/max(sigmaPower);
                gammaPower = powerOutput.avg_Mag(:,6);
                gammaPower = gammaPower/max(gammaPower);    
                
                %xxx move the downsampling up or downsample all the partitioned events
                %
                % downsample the spikeArray
                spikeArray_above0_downsampled = zeros(length(deltaPower), 1);
                spikeArray_1to5_downsampled = zeros(length(deltaPower), 1);
                spikeArray_0to1_downsampled = zeros(length(deltaPower), 1);
                spikeArray_5toMax_downsampled = zeros(length(deltaPower), 1);
                spikeArray_0toneg1_downsampled = zeros(length(deltaPower), 1);
                spikeArray_neg5toMin_downsampled = zeros(length(deltaPower), 1);
                spikeArray_neg1toneg5_downsampled = zeros(length(deltaPower), 1);

                startpoint = 1;
                for ii = 1:length(spikeArray_above0_downsampled)
                    endpoint = startpoint + 1023;
                    spikeArray_above0_downsampled(ii) = sum(spikeArray_above0(startpoint:endpoint));
                    spikeArray_1to5_downsampled(ii) = sum(spikeArray_1to5(startpoint:endpoint));
                    spikeArray_0to1_downsampled(ii) = sum(spikeArray_0to1(startpoint:endpoint));
                    spikeArray_5toMax_downsampled(ii) = sum(spikeArray_5toMax(startpoint:endpoint));
                    spikeArray_0toneg1_downsampled(ii) = sum(spikeArray_0toneg1(startpoint:endpoint));
                    spikeArray_neg1toneg5_downsampled(ii) = sum(spikeArray_neg1toneg5(startpoint:endpoint));
                    spikeArray_neg5toMin_downsampled(ii) = sum(spikeArray_neg5toMin(startpoint:endpoint));
                    startpoint = startpoint + 1024;
                end
                
                % make a random variable for comparison
                    spikeArray_random_downsampled = zeros(length(spikeArray_above0_downsampled), 1);
                    randInd = randsample(1:length(spikeArray_random_downsampled), 200);
                    spikeArray_random_downsampled(randInd) = 1;
               
                % day night breaks
                lightsOn = 1:hoursTicks(11) / 4;
                lightsOff = (hoursTicks(11) / 4):(hoursTicks(23) / 4);
                               
                % plot daytime

                subplot(9, 4, [21, 25])
                
                    lags = [-25:25];
                    powerBand = deltaPower;
                    lightCondition = lightsOn;
                    yylabel = 'Delta Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);

                subplot(9, 4, [22, 26])
                    
                    lags = [-10:10];
                    powerBand = thetaPower;
                    lightCondition = lightsOn;
                    yylabel = 'Theta Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);                     
                    
                    
                subplot(9, 4, [29, 33])
                    
                    lags = [-10:10];
                    powerBand = sigmaPower;
                    lightCondition = lightsOn;
                    yylabel = 'Sigma Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);
                    
                subplot(9, 4, [30, 34])
                    
                    lags = [-25:25];
                    powerBand = gammaPower;
                    lightCondition = lightsOn;
                    yylabel = 'Gamma Power';
                    xxlabel = 'Time (s)';
                    yylim =  [.1, .3];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);

                % plot lightsoff
                subplot(9, 4, [23, 27])
                
                    lags = [-25:25];
                    powerBand = deltaPower;
                    lightCondition = lightsOff;
                    yylabel = 'Delta Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);
                   

                subplot(9, 4, [24, 28])
                     
                    lags = [-10:10];
                    powerBand = thetaPower;
                    lightCondition = lightsOff;
                    yylabel = 'Theta Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);
                   
                subplot(9, 4, [31, 35])
                    
                    lags = [-10:10];
                    powerBand = sigmaPower;
                    lightCondition = lightsOff;
                    yylabel = 'Sigma Power';
                    xxlabel = 'Time (s)';
                    yylim =  [0, .4];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);
                   
                    
                subplot(9, 4, [32, 36])
                
                    lags = [-25:25];
                    powerBand = gammaPower;
                    lightCondition = lightsOff;
                    yylabel = 'Gamma Power';
                    xxlabel = 'Time (s)';
                    yylim =  [.1, .3];
                    
                    patch([-100, 100, 100, -100], [0, 0, 1, 1], 'k', 'FaceColor', 'k', 'FaceAlpha', .10, 'EdgeColor', 'none');
                    hold on
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_random_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_random_downsampled));
                    plot(lags* 4, STA, 'Color', 'k')
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,'k', 'edgecolor', 'k', 'facealpha', .1, 'edgealpha', .1); 
                    
                    xlim([lags(1)*4, lags(end)*4]);
                    ylim(yylim)
                    ylabel(yylabel)
                    xlabel(xxlabel) 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg5toMin_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg5toMin_downsampled));
                    plot(lags* 4, STA, 'Color', colors(1,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(1,:), 'edgecolor', colors(1,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_neg1toneg5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_neg1toneg5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(16,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(16,:), 'edgecolor', colors(16,:), 'facealpha', .1, 'edgealpha', .1);    
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0toneg1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0toneg1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(26,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(26,:), 'edgecolor', colors(26,:), 'facealpha', .1, 'edgealpha', .1);            
                           
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_0to1_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_0to1_downsampled));
                    plot(lags* 4, STA, 'Color', colors(38,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(38,:), 'edgecolor', colors(38,:), 'facealpha', .1, 'edgealpha', .1); 
                    
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_1to5_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_1to5_downsampled));
                    plot(lags* 4, STA, 'Color', colors(48,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(48,:), 'edgecolor', colors(48,:), 'facealpha', .1, 'edgealpha', .1); 
   
                    [STA, STE, ~] = spiketriggeredaverage(lightCondition, spikeArray_5toMax_downsampled(lightCondition), powerBand(lightCondition), lags);
                    STS = std(STE') / sqrt(sum(spikeArray_5toMax_downsampled));
                    plot(lags* 4, STA, 'Color', colors(end,:))
                    x=lags* 4;                      %#initialize x array
                    y1=STA + STS;                   %#create first curve
                    y2=STA - STS;                   %#create second curve
                    X=[x,fliplr(x)];                %#create continuous x value array for plotting
                    Y=[y1,fliplr(y2)];              %#create y values for out and then back
                    fill(X,Y,colors(end,:), 'edgecolor', colors(end,:), 'facealpha', .1, 'edgealpha', .1);
                   
                    %
            % SAVE THE PLOT to its own file
                set(gcf, 'PaperSize', [20 20]) 
                print(gcf, '~/Desktop/Fig 5.pdf', '-dpdf', '-r300')            % clear some things to reset for next plot
            clear SWDsecond SWDStartLOCS keepInd SWDs everySecond SWDInd
                        
        end    
    end
end

%%

% Figure 6

% find example events in wild type animals

% animal 4 is an RR

for i = 4%1:length(F.itemparams)
    for l = 3%1:length(F.itemparams{1, i}.recordDay)
        j = str2num(F.itemparams{1, i}.recordDay{1, l});
        
        % ----- PREPPING DATA FOR PLOTS
        
        % set up some variables from the namefile
        animalID = F.itemparams{1, i}.animalID{1};
        edffilespec = animal(i).day(j).edffilespec;
        [path, name, ext] = fileparts(edffilespec);
        
        % edf for spike triggered averages
        alreadyImportedEDF = read_EDF_mj(edffilespec);
        edf = alreadyImportedEDF.D.edf.signalMat(:, 2);
        fs = alreadyImportedEDF.D.edf.fs;
        
        % time variables for use in plotting
        hourInSeconds = 1:3600;
        totalHours = 1:24;
        twentyfourhours = 86400;
        plotPositions = totalHours - .5;
        
        % labely stuff for use in plotting
        hoursLabels = {'8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '1', '2', '3', '4', '5', '6', '7'};
        hoursTicks = hourInSeconds(end) * (1:24);
        
        % does sleep score and a analyzed file for swds?
        if ~strcmp(F.itemparams{1, i}.sleepScore{1, l}, 'NA') && str2num(F.itemparams{1, i}.SWDdetectOK{1, l}) == 1
        
            % set up sleep stage stuff
            sleepScoresText = animal(i).day(j).sleepScore.score;
            wakeInd = find(strcmp(sleepScoresText, 'Wake'));
            REMInd = find(strcmp(sleepScoresText, 'REM'));
            nREMInd = find(strcmp(sleepScoresText, 'Non REM'));
            artInd = find(strcmp(sleepScoresText, 'Artifact'));
            sleepScores = zeros(length(sleepScoresText), 1);
            sleepScores(wakeInd) = 1;
            sleepScores(REMInd) = 3;
            sleepScores(nREMInd) = 2;
            sleepScores(artInd) = NaN;
            time = animal(i).day(j).sleepScore.time;         
            
            % setup up some stuff for plotting SWDs
            SWDs = animal(i).day(j).SWDs;

            for ii = 1:size(SWDs, 2)
                SWDStartLOCS(ii) = SWDs(ii).SWDLOCS(1); 
                
            end
           
            %{
            % filter out all the ones based on the humanScored index that were 100 percent not SWDs as manually decided.
            SWDStartLOCS = SWDStartLOCS(find(animal(i).day(j).classificationResponse == 1));
            nonSWDStartLOCS = SWDStartLOCS(find(animal(i).day(j).classificationResponse == 2));
            %}
            
            % now expand this record to a mark for every second in which an SWD was triggered.
            SWDStartLOCS = floor(SWDStartLOCS/256); % convert to seconds instead of sample points
            everySecond = 1:time(end); % make a vector that is 24 hours worth in seconds
            SWDsecond = zeros(1, size(everySecond, 2))'; % initialize ones for SWD presence at these seconds
            for ii = 1:length(SWDStartLOCS)
                SWDsecond(SWDStartLOCS(ii)) = 1; % now add a SWD where they occur in seconds
            end
            
            % drop pairs of SWDs that start in the same second..
            if sum(SWDsecond) ~= length(SWDStartLOCS)
               diffSWDStartLOCS = [ 10, diff(SWDStartLOCS)]; % 10 is just a placeholder bigger than the next threshold
               keepInd = find(diffSWDStartLOCS);
               SWDStartLOCS = SWDStartLOCS(keepInd);
            else
                keepInd = 1:sum(SWDsecond);
            end
               
            
            % expand the sleep scores from epochs to an every second record. everySecondWake, Sleep and nREM
            sleepScoresExpanded = [];
            for ii = 1:length(sleepScores) - 1
                sleepScoresExpanded((ii * 4)) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 1) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 2) = sleepScores(ii);
                sleepScoresExpanded((ii * 4) - 3) = sleepScores(ii);
            end
               
            % find the intersect of SWDs occuring at each sleep score and produce an index
            wiswd = sleepScoresExpanded == 1 & SWDsecond' == 1;
            riswd = sleepScoresExpanded == 3 & SWDsecond' == 1;
            nriswd = sleepScoresExpanded == 2 & SWDsecond' == 1;
            
            % make a matrix of the SWD characteristics to plot intensity -- perhaps amplitude and duration?
            SWDInd = find(SWDsecond);
            critsArray = zeros(1, size(everySecond, 2))'; % this is a zeros index so that 
            classificationScore = animal(i).day(j).classificationScore;
            classificationScore = classificationScore(keepInd);
            classificationResponse = animal(i).day(j).classificationResponse;
            classificationResponse = classificationResponse(keepInd);
            
            
            % add negative value for events that are NOT SWDs
            for ii = 1:length(keepInd)
                if classificationResponse(ii) == 1
                    classificationScore(ii) = -classificationScore(ii);
                end
            end
            
            critsArray(SWDInd) = classificationScore;

           
                
            
            % SWD AVERAGES PER HOUR
            yesYes = SWDsecond & (critsArray > 0); % you can change this threshold to a more negative number to include more events and a more positive number to include fewer events.
            
            WAKE_SWDsecond = yesYes' == 1 & wiswd == 1;
            REM_SWDsecond = yesYes' == 1 & riswd == 1;
            nREM_SWDsecond = yesYes' == 1 & nriswd == 1;
            for ii = totalHours
                hourInd = hourInSeconds + (hourInSeconds(end) * ii) - hourInSeconds(end);
                seizuresPerHour(ii) = sum(yesYes(hourInd));
                WAKE_seizuresPerHour(ii) = sum(WAKE_SWDsecond(hourInd));
                REM_seizuresPerHour(ii) = sum(REM_SWDsecond(hourInd));
                nREM_seizuresPerHour(ii) = sum(nREM_SWDsecond(hourInd));
            end
            
            % SLEEP STAGE AVERAGES PER HOUR
            for ii = totalHours
                hourInd = hourInSeconds + (hourInSeconds(end) * ii) - hourInSeconds(end);
                WAKE_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 1);
                percentage_WAKE_perHour(ii) = WAKE_perHour(ii);
                REM_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 3);
                percentage_REM_perHour(ii) = REM_perHour(ii);
                nREM_perHour(ii) = sum(sleepScoresExpanded(hourInd) == 2);
                percentage_nREM_perHour(ii) = nREM_perHour(ii); 
            end
            % convert to percentages
            percentage_WAKE_perHour = percentage_WAKE_perHour / hourInSeconds(end);
            percentage_REM_perHour = percentage_REM_perHour / hourInSeconds(end);
            percentage_nREM_perHour = percentage_nREM_perHour / hourInSeconds(end);


        end
    end
end

%% 
% use the above code to pull the locations of the seizures

seizureIntensity = critsArray(find(critsArray > 0));
seizureLocs = find(critsArray > 0) * 256;

signal = normalizeEEG(edf, 256);

%% figure for animal 4
figure('units', 'inch', 'pos', [10 10 4 8]);
subplot(6, 1, 1)
    plot(signal(seizureLocs(10)-300:seizureLocs(10)+ 700), 'k')
    axis off
    ylim([-10, 10])
subplot(6, 1, 2)
    plot(signal(seizureLocs(21)-400:seizureLocs(21)+ 600), 'k')
    axis off
    ylim([-10, 10])
subplot(6, 1, 3)
    plot(signal(seizureLocs(27)-500:seizureLocs(27)+ 500), 'k')
    axis off
    ylim([-10, 10])
subplot(6, 1, 4)
    plot(signal(seizureLocs(8)-300:seizureLocs(8)+ 700), 'k')
    axis off
    ylim([-10, 10])
subplot(6, 1, 5)
    plot(signal(seizureLocs(5)-400:seizureLocs(5)+ 600), 'k')
    axis off
    ylim([-10, 10])
subplot(6, 1, 6)
    plot(signal(seizureLocs(47)-300:seizureLocs(47)+ 700), 'k')
    axis off
    ylim([-10, 10])

                    set(gcf, 'PaperSize', [20 20]) 
                print(gcf, '~/Desktop/Fig 6.pdf', '-dpdf', '-r300')   
