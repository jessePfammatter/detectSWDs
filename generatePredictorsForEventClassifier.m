function [ featureSet ] = generatePredictorsForEventClassifier( eventsObject, fs )
% generatePredictorsForEventClassifier calculates 18 characteristics from
% and EEG signal. This script takes an 'eventsObject' as generated in the 
% 'detectSWDs_automated script.m' script. The characteristics are the event
% duration (length of the signal, in sample points), line length, and
% normalized line length of the signal along with the mean, std, max, line 
% length and, normalized line length of 3 vectors of data complied from 
% wavelet transforms. The three vectors are designed to target the 
% freqeuncy characteristics of Spike-and-Wave discharges and are as 
% follows:
%
% sixHzSignal: the sum of the wavelet (Mortlet) decomposition signals aggregated 
% accross scales 40:48. This is designed to hone in on the ~6Hz
% characteristic of the spike-and-wave discharges (SWDs).
%
% harmonicSignal: the sum of the wavelet (Mortlet) decomposition signals aggregated 
% accross scales 20:28. This is designed to hone in on the ~16 - 32Hz
% plusing harmonic signal that we see in SWDs.
%
% DB4: the sum of the wavelet decopmosition (DB4) aggregated accross the 
% scales 25:45. This is designed to hone in on the ~6Hz characteristics of
% the SWDs but using the DB4 wavelet which has been suggested to be useful
% in identifying 'spikes' given it's sharp features.
%
% For an example, use 'detectSWDs_automated.m' to create and eventsObject
% in order to test this script.
%
% JP 2017

    % length of the singal clip used
    lengthClip = fs * 8; % for an 8 second time window
    
    % a buffer around the exact start and stop time (in sample points) that seems relevant to the SWD.
    buffer = floor(fs/25);

    % calculate start and stop times of events -- these are not perfect but will be correct later in the script with spike detection stuff
    for i = 1:size(eventsObject, 2)
        eventStart(i) = floor((.5*lengthClip) - (.5*eventsObject(i).seizureDuration)) - buffer;
        if eventStart(i) < 1 + buffer
            eventStart(i) = 1 + buffer;
        end
        eventStop(i) = ceil((.5*lengthClip) + (.5*eventsObject(i).seizureDuration)) + buffer;
        if eventStop(i) > lengthClip - buffer
            eventStop(i) = lengthClip - buffer;
        end
    end

    for i = 1:size(eventsObject, 2)    

        % wavelet analyses on the full signal clips
        eventsObject(i).cwtMortOut = cwt(eventsObject(i).signalClips, 'amor', fs);

        % sixHz variables
        eventsObject(i).sixHz = sum(abs(eventsObject(i).cwtMortOut(floor(fs/6.4):floor(fs/5.22),:))); %5-8hz approximately
        sixHz_mean(i) = mean(eventsObject(i).sixHz(eventStart(i):eventStop(i)));
        sixHz_std(i) = std(eventsObject(i).sixHz(eventStart(i):eventStop(i)));
        sixHz_max(i) = max(eventsObject(i).sixHz(eventStart(i):eventStop(i)));

        % inBetween variables
        eventsObject(i).inBetween = sum(abs(eventsObject(i).cwtMortOut(floor(fs/8.53):floor(fs/6.56),:))); 
        inBetween_mean(i) = mean(eventsObject(i).inBetween(eventStart(i):eventStop(i)));
        inBetween_std(i) = std(eventsObject(i).inBetween(eventStart(i):eventStop(i)));
        inBetween_max(i) = max(eventsObject(i).inBetween(eventStart(i):eventStop(i)));
        
        % calculate the harmonic  variables
        eventsObject(i).harmonic = sum(abs(eventsObject(i).cwtMortOut(floor(fs/12.8):floor(fs/8.83),:))); %15-32 hrz approx
        harmonic_mean(i) = mean(eventsObject(i).harmonic(eventStart(i):eventStop(i)));
        harmonic_std(i) = std(eventsObject(i).harmonic(eventStart(i):eventStop(i)));
        harmonic_max(i) = max(eventsObject(i).harmonic(eventStart(i):eventStop(i)));
        
        % higherFreq variables
        eventsObject(i).higherFreq = sum(abs(eventsObject(i).cwtMortOut(floor(fs/25.6):floor(fs/13.47),:))); 
        higherFreq_mean(i) = mean(eventsObject(i).higherFreq(eventStart(i):eventStop(i)));
        higherFreq_std(i) = std(eventsObject(i).higherFreq(eventStart(i):eventStop(i)));
        higherFreq_max(i) = max(eventsObject(i).higherFreq(eventStart(i):eventStop(i)));


    end

    % redistribute to the featureSet object
    featureSet = table(sixHz_mean', sixHz_std', sixHz_max', inBetween_mean', inBetween_std', inBetween_max', harmonic_mean',  harmonic_std', harmonic_max', higherFreq_mean', higherFreq_std', higherFreq_max');
    featureSet.Properties.VariableNames = {'sixHz_mean', 'sixHz_std', 'sixHz_max', 'inBetween_mean', 'inBetween_std', 'inBetween_max', 'harmonic_mean', 'harmonic_std', 'harmonic_max',  'higherFreq_mean', 'higherFreq_std', 'higherFreq_max'};

end
