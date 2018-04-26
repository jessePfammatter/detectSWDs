function varargout = eeggui(varargin)
% EEGGUI MATLAB code for eeggui.fig
%      EEGGUI, by itself, creates a new EEGGUI or raises the existing
%      singleton*.
%
%      H = EEGGUI returns the handle to a new EEGGUI or the handle to
%      the existing singleton*.
%
%      EEGGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EEGGUI.M with the given input arguments.
%
%      EEGGUI('Property','Value',...) creates a new EEGGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before eeggui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eeggui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eeggui

% Last Modified by GUIDE v2.5 16-Oct-2016 20:22:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @eeggui_OpeningFcn, ...
                   'gui_OutputFcn',  @eeggui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before eeggui is made visible.
function eeggui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eeggui (see VARARGIN)

% Choose default command line output for eeggui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% ******** USER DEFAULTS **********
defaultfig
axes(handles.EDF_Axes)
    xlabel('Time, s')
    ylabel('Signals, arbitrary units')

handles.D.edf.npts = 32 .* 256;  % sample pts
handles.D.edf.currpt = 1;       % current sample pt

handles.D.edf.prepostprocdata.dumpcount = 0;
handles.D.edf.prepostprocdata.dump = {}; % {timept, delta power, motion},etc - temp analysis storage

set(handles.Start_Time_Edit, 'str', '--')
set(handles.Length_Edit, 'str', '--')


% ******** USER DEFAULTS **********


% Update handles structure
guidata(hObject, handles);
assignin('base','handles',handles)
% UIWAIT makes eeggui wait for user response (see UIRESUME)
% uiwait(handles.EEG_GUI);


% --- Outputs from this function are returned to the command line.
function varargout = eeggui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Open_EDF_Callback(hObject, eventdata, handles)
% hObject    handle to Open_EDF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [f p] = uigetfile('*.edf', 'Open EDF file:');
    cd(p)
    handles.D.edf.filespec = [p f];
%     disp(['Please wait while opening EDF file: ' handles.D.edf.filespec])
    
%     % ********* The old way ***********
%     [handles.D.edf.header handles.D.edf.signalHeader signalCell] = blockEdfLoad(handles.D.edf.filespec);
% %     EDFdummy = read_EDF_mj(handles.D.edf.filespec);
%     handles.D.edf.signalMat = cell2mat(signalCell);
%     handles.D.edf.signalMax = max( handles.D.edf.signalMat );
%     handles.D.edf.signalMin = min( handles.D.edf.signalMat );
%     handles.D.edf.signalStd = std( handles.D.edf.signalMat );
%     handles.D.edf.nchans = handles.D.edf.header.num_signals;
%     handles.D.edf.nrecs = handles.D.edf.header.num_data_records;
%     handles.D.edf.recdur = handles.D.edf.header.data_record_duration;
%     handles.D.edf.nsecs = handles.D.edf.nrecs .* handles.D.edf.recdur;
%     dummy = cat(1, handles.D.edf.signalHeader );
%     dummy = cat(1, dummy.samples_in_record);
%     if length(unique(dummy)) == 1
%         handles.D.edf.fs = unique(dummy);   % sampling rate, Hz
%     else
%         warndlg({'NOTE: There are different sampling rates on different channels.'; 'I''m using the first one.'; 'This may cause errors later.'}, 'Sampling Mismatch!')
%         handles.D.edf.fs = dummy(1);        % sampling rate, Hz (first channel only)
%     end  
%     headertext = {  handles.D.edf.filespec;...
%                     ['sampling rate = ' num2str(handles.D.edf.fs) ' Hz'];...
%                     ['total time = ' num2str(handles.D.edf.nrecs .* handles.D.edf.recdur) ' seconds  (' num2str([handles.D.edf.nrecs .* handles.D.edf.recdur ./ 60 ./ 60], '%2.2f') ' hours) '] ...
%                   };   
%     set(handles.Header_Info, 'str', headertext)
    % ********* The old way ***********
    
    % ********* The new way, that uses read_EDF_mj.m ***********
    %   Some code duplicates that in read_EDF_mj to avoid overwriting the
    %   handles variable.
    EDFdummy = read_EDF_mj(handles.D.edf.filespec);
    handles.D.edf.eegchans  = EDFdummy.D.edf.eegchans;   % Channels with eeg signals
    handles.D.edf.epochdur  = EDFdummy.D.edf.epochdur;     % sec
    handles.D.edf.currpt    = EDFdummy.D.edf.currpt;       % current sample pt
    handles.D.edf.signalMat = EDFdummy.D.edf.signalMat;
    handles.D.edf.eegData   = EDFdummy.D.edf.eegData;       % This is a possible duplicate of what is in signalMat. If so, should fix to not duplicate huge data. But watch out for breaking other code elsewhere
    handles.D.edf.signalMax = max( EDFdummy.D.edf.signalMat );
    handles.D.edf.signalMin = min( EDFdummy.D.edf.signalMat );
    handles.D.edf.signalStd = std( EDFdummy.D.edf.signalMat );
    handles.D.edf.header        = EDFdummy.D.edf.header;
    handles.D.edf.signalHeader  = EDFdummy.D.edf.signalHeader;
    handles.D.edf.num_signals           = EDFdummy.D.edf.header.num_signals;
    handles.D.edf.num_data_records      = EDFdummy.D.edf.header.num_data_records;
    handles.D.edf.data_record_duration 	= EDFdummy.D.edf.header.data_record_duration;
    handles.D.edf.nchans    = EDFdummy.D.edf.header.num_signals;
    handles.D.edf.nrecs     = EDFdummy.D.edf.header.num_data_records;
    handles.D.edf.recdur    = EDFdummy.D.edf.header.data_record_duration;
    handles.D.edf.nsecs     = EDFdummy.D.edf.nrecs .* handles.D.edf.recdur;
    dummy = cat(1, EDFdummy.D.edf.signalHeader );
    dummy = cat(1, dummy.samples_in_record);
    if length(unique(dummy)) == 1
        handles.D.edf.fs = unique(dummy);   % sampling rate, Hz
    else
        warndlg({'NOTE: There are different sampling rates on different channels.'; 'I''m using the first one.'; 'This may cause errors later.'}, 'Sampling Mismatch!')
        handles.D.edf.fs = dummy(1);        % sampling rate, Hz (first channel only)
    end  
    headertext = {  handles.D.edf.filespec;...
                    ['sampling rate = ' num2str(EDFdummy.D.edf.fs) ' Hz'];...
                    ['total time = ' num2str(handles.D.edf.nrecs .* handles.D.edf.recdur) ' seconds  (' num2str([handles.D.edf.nrecs .* handles.D.edf.recdur ./ 60 ./ 60], '%2.2f') ' hours) '] ...
                  };   
    set(handles.Header_Info, 'str', headertext)
    % ********* The new way, that uses read_EDF_mj.m ***********
    
    guidata(hObject, handles);
    assignin('base','handles',handles)
    Update_Callback(hObject, eventdata, handles)
    disp('  DONE')
    chirpsound_mj

    
function Draw_EDF_Data(hObject, eventdata, handles)
    nchans = handles.D.edf.nchans;       % sample pts
    npts = handles.D.edf.npts;       % sample pts
    currpt = handles.D.edf.currpt;       % current sample pt
    lastpt = size(handles.D.edf.signalMat, 1);

    fs = handles.D.edf.fs;
    fnyquist = fs/2;                        %  Nyquist frequency
    bin_vals = [0 : npts-1];              	% Freq bins
    fax_Hz = bin_vals*fs/npts;               % Freq Hz
    N_2 = ceil(npts/2);                     % Single-sided PSD 
    
    bands.delta.range = [0.5 4];                  % Freq ranges for bands to analyze
    bands.SWD.range   = [4.1 6.9];                    
    bands.theta.range = [7.1 12];   
    bands.gamma.range = [25 100];
    
    time = [currpt:currpt+npts-1] ./ fs;          % seconds
    tracescale = 0.1;                              
%     disp([currpt npts])
    handles.D.edf.currdata = handles.D.edf.signalMat(currpt:currpt+npts-1, :);    % clip out appropriate section of data
    dispdata = handles.D.edf.currdata;
    
    
%   Pre/Post Processing
    if isfield(handles.D.edf, 'prepostprocsettings')
        handles.D.edf.prepostprocdata.Hz  = fax_Hz(1:N_2)';     % Freq Axis     
        handles.D.edf.prepostprocdata.PSD =[];
        for chan = 1:length(handles.D.edf.prepostprocsettings.filterflag);
%             Prefilter data to DISPLAY 
            if handles.D.edf.prepostprocsettings.filterflag(chan)
                dispdata(:, chan) = filter(handles.D.edf.prepostprocsettings.Hd(chan), dispdata(:, chan));
            end
%             Compute PSDs on UNFILTERED data
            if handles.D.edf.prepostprocsettings.PSDflag(chan)
%                 X_mags = abs(fft(handles.D.edf.currdata(:,chan)));          % Power Spectrum
%                 handles.D.edf.prepostprocdata.PSD(:,chan) = X_mags(1:N_2);
                [X_mags,f] = pwelch(handles.D.edf.currdata(:,chan), 960, 100, logspace(log10(0.5), log10(100), 64), fs);
                handles.D.edf.prepostprocdata.Hz = f(:);
                handles.D.edf.prepostprocdata.PSD(:,chan) = X_mags(:);
                
                bandnames = fieldnames(bands);
                for band = 1:length(bandnames)
                    rangelims = eval(['bands.' bandnames{band} '.range' ]);
                    indxs = find( f>=rangelims(1) &  f<=rangelims(2) );
%                     avg_freq(band)    = mean( f( indxs ) ) ;  
%                     avg_Mag(chan, band)     = mean( X_mags( indxs ) );
                        [max_Mag(band), mxindx]     = max( X_mags( indxs ) );
                        max_freq(band)              = f(indxs(mxindx)) ;
                        avg_Mag(band)               = mean( X_mags( indxs ) );
                        avg_freq(band)              = mean(f(indxs)) ;                
                end                   
                handles.D.edf.prepostprocdata.channel(chan).max_freq    = max_freq;
                handles.D.edf.prepostprocdata.channel(chan).max_Mag     = max_Mag;
                handles.D.edf.prepostprocdata.channel(chan).avg_freq    = avg_freq;
                handles.D.edf.prepostprocdata.channel(chan).avg_Mag     = avg_Mag;  
            end
            

        end    
    end

%     Scale EEG traces for display
    datarange = handles.D.edf.signalMax - handles.D.edf.signalMin;
    dispdata = dispdata ./ repmat(handles.D.edf.signalStd, size(dispdata, 1), 1);
    dispdata = dispdata .* tracescale + repmat( [1:nchans], size(dispdata, 1), 1); 

%     Draw EEG Traces
    axes(handles.EDF_Axes)
    plot(time, dispdata )
    set(handles.EDF_Axes, 'ydir', 'rev')
    set(handles.EDF_Axes, 'xlim', [currpt currpt+npts]./fs, 'ylim', [0 nchans+1], 'ytick', [1:nchans], 'yticklabel', num2str([1:nchans]'))
    xlabel('Time, s')
    ylabel('Channel #')     
    
%     Draw PSDs    
    if isfield(handles.D.edf, 'prepostprocsettings') & any( handles.D.edf.prepostprocsettings.PSDflag)
        axes(handles.FFT_Segment)
        cla
        plot(handles.D.edf.prepostprocdata.Hz, handles.D.edf.prepostprocdata.PSD )
        hold on
        xlabel('Frequency (Hz)')
        ylabel('Magnitude');
        set(gca, 'xsc', 'log', 'ysc', 'lin')
        set(gca, 'xlim', [handles.D.edf.prepostprocdata.Hz([1 end])])
%       
        colord = get( handles.FFT_Segment, 'ColorOrder');
        for chan = 1:size(handles.D.edf.prepostprocdata.PSD, 2)
            plot(handles.D.edf.prepostprocdata.channel(chan).max_freq, handles.D.edf.prepostprocdata.channel(chan).max_Mag, 'o', 'color', colord(chan, :))
            plot(handles.D.edf.prepostprocdata.channel(chan).avg_freq, handles.D.edf.prepostprocdata.channel(chan).avg_Mag, 'o', 'color', colord(chan, :), 'markerfacecolor', colord(chan, :))
        end
        

% START HERE 
% - FIGURE OUT HOW TO KEEP TRACK OF CHANNEL & BAND SCORES            
% - Take averges across CHANNELS, and plot, and save.        
% - Plot MOTION versus power in various bands - Maybe PCA??
        dummy = cat(1, handles.D.edf.prepostprocdata.channel);
        x = mean( cat(1, dummy.avg_freq)) ;
        y = mean( cat(1, dummy.avg_Mag)) ;
        x2 = mean( cat(1, dummy.max_freq)) ;
        y2 = mean( cat(1, dummy.max_Mag)) ;
        plot(x, y, 'kv-', 'markerfacecolor', 'k', 'markersize', 12)
        plot(x2, y2, 'r^-', 'markerfacecolor', 'r', 'markersize', 12)
        
        lineord = {'-', '--', '-.', ':'};
        ylim = get(gca, 'ylim');
        for band = 1:length(bandnames)      
            rangelims = eval(['bands.' bandnames{band} '.range' ]);
            plot([rangelims; rangelims], [ylim', ylim'], 'k', 'linestyle', lineord{band})
        end   
        hold off
        
        
    % MOTION        
        motionchan = find(handles.D.edf.prepostprocsettings.motionflag);
    	locomotion = mean(handles.D.edf.currdata(:,motionchan));

        
        handles.D.edf.prepostprocdata.dumpcount = handles.D.edf.prepostprocdata.dumpcount+1;
        epoch = handles.D.edf.prepostprocdata.dumpcount;
        for chan = 1:size(handles.D.edf.prepostprocdata.PSD, 2)
            for band = 1:length(bandnames)
                
                handles.D.edf.prepostprocdata.dump.bandnames = bandnames;
                handles.D.edf.prepostprocdata.dump.EDFstats.max_freq{chan, band, epoch}     = handles.D.edf.prepostprocdata.channel(chan).max_freq(band);
                handles.D.edf.prepostprocdata.dump.EDFstats.max_Mag{chan, band, epoch}      = handles.D.edf.prepostprocdata.channel(chan).max_Mag(band);
                handles.D.edf.prepostprocdata.dump.EDFstats.avg_freq{chan, band, epoch}     = handles.D.edf.prepostprocdata.channel(chan).avg_freq(band);
                handles.D.edf.prepostprocdata.dump.EDFstats.avg_Mag{chan, band, epoch}      = handles.D.edf.prepostprocdata.channel(chan).avg_Mag(band);
                handles.D.edf.prepostprocdata.dump.EDFstats.locomotion(epoch)               = locomotion;
                handles.D.edf.prepostprocdata.dump.EDFstats.starttime(epoch)                = currpt./fs;
                handles.D.edf.prepostprocdata.dump.EDFstats.epochdur(epoch)                 = npts./fs;
            
            end
        end   
        
    end
    
    guidata(hObject, handles);
    assignin('base','handles', handles) 
%     disp(['Just ran Draw_EDF_Data ' num2str(handles.D.edf.prepostprocdata.dumpcount)])
    
    
    
    

% --- Executes on button press in Update.
function Update_Callback(hObject, eventdata, handles)
    set(handles.Start_Time_Edit, 'str', num2str(handles.D.edf.currpt./handles.D.edf.fs, '%2.4f') )
    set(handles.Length_Edit, 'str', num2str(handles.D.edf.npts./handles.D.edf.fs) )
    Draw_EDF_Data(hObject, eventdata, handles)
    handles = guidata(hObject);
    guidata(hObject, handles);
    assignin('base','handles',handles) 
    disp(['Updated ' num2str(handles.D.edf.prepostprocdata.dumpcount)])

% --- Executes on key press with focus on EEG_GUI or any of its controls.
function EEG_GUI_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to EEG_GUI (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
    npts = handles.D.edf.npts;       % sample pts
    currpt = handles.D.edf.currpt;       % current sample pt
    lastpt = size(handles.D.edf.signalMat, 1);
%     disp( [currpt npts currpt+npts lastpt] )
    switch eventdata.Key
        case 'rightarrow'   % Move to next time window
            if currpt+npts >= lastpt;
                warndlg('End of Data', 'Warning', 'modal')
            else
                currpt = ceil( currpt + npts );
            end    
        case 'leftarrow'    % Move to previous time window
            if currpt-npts < 1;
                warndlg('Beginning of Data', 'Warning', 'modal')
            else
                currpt = ceil( currpt - npts );
            end
        case 'uparrow'      % Double time window
            npts = floor( 2 .* npts );
        case 'downarrow'    % Halve time window
            npts = floor( 0.5 .* npts );
        case 'comma'    	% Go to start
            currpt = 1;
        case 'period'        % Go to end
            currpt = ceil( lastpt - npts );
    end       
    handles.D.edf.npts = npts;       % sample pts
    handles.D.edf.currpt = currpt;       % current sample pt
    
    Update_Callback(hObject, eventdata, handles)
    handles = guidata(hObject);

    guidata(hObject, handles);
    assignin('base', 'handles', handles);        



function Start_Time_Edit_Callback(hObject, eventdata, handles)    
    npts = handles.D.edf.npts;       % sample pts
    currpt = handles.D.edf.currpt;       % current sample pt
    lastpt = size(handles.D.edf.signalMat, 1);
    fs = handles.D.edf.fs;
    testpt =  str2num( get(hObject, 'str') ) .* fs;
    if testpt < 1; testpt = 1; end
    if testpt+npts > lastpt;
                warndlg('End of Data', 'Warning', 'modal')
    elseif testpt < 1
                warndlg('Beginning of Data', 'Warning', 'modal')
    else
        handles.D.edf.currpt = ceil( testpt );
    end
    Update_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
    assignin('base', 'handles', handles);        
    
% --- Executes during object creation, after setting all properties.
function Start_Time_Edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Length_Edit_Callback(hObject, eventdata, handles)
    npts = handles.D.edf.npts;       % sample pts
    currpt = handles.D.edf.currpt;       % current sample pt
    lastpt = size(handles.D.edf.signalMat, 1);
    fs = handles.D.edf.fs;
    testpt =  str2num( get(hObject, 'str') ) .* fs
    if currpt+testpt >= lastpt;
                warndlg('End of Data', 'Warning', 'modal')
    elseif testpt < 1
                warndlg('Bad Value', 'Warning', 'modal')
    else
        handles.D.edf.npts = testpt;
    end
    Update_Callback(hObject, eventdata, handles)
    guidata(hObject, handles);
    assignin('base', 'handles', handles);        
    
% --- Executes during object creation, after setting all properties.
function Length_Edit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Help.
function Help_Callback(hObject, eventdata, handles)
    helpstr = {
        'MENUS:';...
        'File Menu      -  Open, Save, etc';...
        '   ';...
        'EDIT BOXES';...
        'Time     -  Start of data window to display, in seconds.'; ...
        'Dur      -  Duration of data window to display, in seconds.'; ...
        '   '   ;...
        'BUTTONS (These are a hangover from a previous GUI -'; ...
        '            and need to be customized)';...
        'Mark           -  Bookmark the current window'; ...
        'Next Mark      -  Jump ahead to the next bookmarked view'; ...
        'Previous Mark  -  Jump back to the previous bookmarked view'; ...
        'Clear Marks    -  Delete all bookmarks'; ...
        'Export Marks   -  Save a MAT file with the current marks, '; ...
        '                       and an AVI file with the current marked views'; ...
        'Import Marks   -  Read a MAT file with previously saved bookmarks'; ...
        ''; ...
        'Keyboard shortcuts:';...
        '                   RIGHTARROW  - Jump ahead';...
        '                   LEFTARROW   - Jump back';...
        '                   UPARROW     - Expand window length';...
        '                   DOWNARROW   - Shrink window length';...
        '                   < (,)       - Jump to start of data';...
        '                   > (.)       - Jump to end of data';...

        }
    out = dialog('WindowStyle', 'normal', 'Name', 'My Dialog');
    axes('units', 'norm', 'pos', [0 0 1 1])
    axis off
    text(0.5, 0.95, 'EEG GUI Help', 'fontname', 'courier', 'fontsize', 18, 'horiz', 'center')
    text(0.05, 0.5, helpstr, 'fontname', 'courier', 'fontsize', 14)


% --------------------------------------------------------------------    
% Cancel Pushbutton callback
function prepostprocsettings_cancel_call(hObject, eventdata, handles)
    delete(handles.D.edf.prepostprocsettings.figh)
    
    
% --------------------------------------------------------------------
function Analysis_Callback(hObject, eventdata, handles)



% --------------------------------------------------------------------
function Pre_Post_Processing_Callback(hObject, eventdata, handles)
    nchans = handles.D.edf.nchans;       % sample pts
    % Create figure
        handles.D.edf.prepostprocsettings.figh = figure('units','inch','position',[3 5 6 6], 'toolbar','none','menu','none', 'name', 'Pre/Post Processing Settings');
        % Create yes/no checkboxes
        xoff1    = 0.01;    
        wid1     = 0.2;
        uicontrol('style','text','units','norm', 'position',[0.01 0.9 0.5 0.1],'string','Apply Pre/Post Processing:', 'fontname', 'times', 'fontsize', 18, 'back', 'w');
        uicontrol('style','text','units','norm', 'position',[xoff1+0.4 0.8 0.1 0.1],'string', 'F(low)', 'fontname', 'times', 'fontsize', 14, 'back', 'w');
        uicontrol('style','text','units','norm', 'position',[xoff1+0.5 0.8 0.1 0.1],'string', 'F(hi)',  'fontname', 'times', 'fontsize', 14, 'back', 'w');
        uicontrol('style','text','units','norm', 'position',[xoff1+0.675 0.8 0.1 0.1],'string', 'PSD',  'fontname', 'times', 'fontsize', 14, 'back', 'w');
        uicontrol('style','text','units','norm', 'position',[xoff1+0.775 0.8 0.1 0.1],'string', 'Motion',  'fontname', 'times', 'fontsize', 14, 'back', 'w');

        for chan = 1:nchans
            ybase = 0.95; yoff = 0.85./nchans; hgt = 0.3./nchans;
            handles.D.edf.prepostprocsettings.channame(chan) = uicontrol('style','text','units','norm', 'position',[xoff1 ybase-(chan.*yoff) wid1 hgt],  'string',   ['Channel #' num2str(chan)],  'value', 1, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
            handles.D.edf.prepostprocsettings.filtercheck(chan) = uicontrol('style','checkbox','units','norm', 'position',[xoff1+0.2 ybase-(chan.*yoff) wid1 hgt],  'string',   ['Prefilter'],  'value', 0, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
            handles.D.edf.prepostprocsettings.filterlow(chan) = uicontrol('style','edit','units','norm', 'position',[xoff1+0.4 ybase-(chan.*yoff) wid1./2 hgt],  'string',   ['0.5'],  'value', 1, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
            handles.D.edf.prepostprocsettings.filterhigh(chan) = uicontrol('style','edit','units','norm', 'position',[xoff1+0.5 ybase-(chan.*yoff) wid1./2 hgt],  'string',   ['100'],  'value', 1, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
            handles.D.edf.prepostprocsettings.PSDcheck(chan) = uicontrol('style','checkbox','units','norm', 'position',[xoff1+0.7 ybase-(chan.*yoff) wid1./2 hgt],  'string',   [''],  'value', 1, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
            handles.D.edf.prepostprocsettings.motioncheck(chan) = uicontrol('style','checkbox','units','norm', 'position',[xoff1+0.8 ybase-(chan.*yoff) wid1./2 hgt],  'string',   [''],  'value', 0, 'fontname', 'times', 'fontsize', 14, 'back', 'w');
        end
        
        
        uicontrol('style','pushbutton', 'units','norm', 'position',[xoff1+0.6 ybase-(chan.*yoff)-0.1 0.8.*wid1 hgt],         'string','Cancel','callback',   {@prepostprocsettings_cancel_call, handles} );
        uicontrol('style','pushbutton', 'units','norm', 'position',[xoff1+0.8 ybase-(chan.*yoff)-0.1 0.8.*wid1 hgt],         'string','OK','callback',       {@prepostprocsettings_ok_call, handles} );
    guidata(handles.EEG_GUI, handles);
    assignin('base', 'handles', handles);        
 
    
function prepostprocsettings_ok_call(hObject, eventdata, handles)
    disp('OK callback')
    vals = get(handles.D.edf.prepostprocsettings.filtercheck, 'Value');     % Check whether to prefilter
        handles.D.edf.prepostprocsettings.filterflag = cat( 1, vals{:} );
    vals = get(handles.D.edf.prepostprocsettings.PSDcheck, 'Value');     % Check whether to compute PSD
        handles.D.edf.prepostprocsettings.PSDflag = cat( 1, vals{:} );
    vals = get(handles.D.edf.prepostprocsettings.motioncheck, 'Value');     % Check whether to compute Motion score
        handles.D.edf.prepostprocsettings.motionflag = cat( 1, vals{:} );        

    %     Create bandpass pre-filters for selected channels
    Fs = handles.D.edf.fs;  % Sampling Freq
    for chan = 1:length(handles.D.edf.prepostprocsettings.filterflag)
        if handles.D.edf.prepostprocsettings.filterflag(chan)
            Fstop1 = 0.95.*str2num(get(handles.D.edf.prepostprocsettings.filterlow(chan), 'string')); % First Stopband Frequency
            Fpass1 = str2num(get(handles.D.edf.prepostprocsettings.filterlow(chan), 'string')); % First Passband Frequency
            Fpass2 = str2num(get(handles.D.edf.prepostprocsettings.filterhigh(chan), 'string')); % Second Passband Frequency
            Fstop2 = 1.05.*str2num(get(handles.D.edf.prepostprocsettings.filterhigh(chan), 'string')); % Second Stopband Frequency
            Astop1 = 80; % First Stopband Attenuation (dB)
            Apass = 1; % Passband Ripple (dB)
            Astop2 = 80; % Second Stopband Attenuation (dB)
            match = 'stopband'; % Band to match exactly
            % Construct an FDESIGN object and call its BUTTER method.
            h = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
            handles.D.edf.prepostprocsettings.Hd(chan) = design(h, 'butter', 'MatchExactly', match);
        end    
    end
    delete(handles.D.edf.prepostprocsettings.figh)
%     Update_Callback(hObject, eventdata, handles)
%     handles = guidata(hObejct);
    guidata(handles.EEG_GUI, handles);
    assignin('base', 'handles', handles);

   

    
    


    
 


% --------------------------------------------------------------------
function Detect_SWDs_Callback(hObject, eventdata, handles)
    EDF.D.edf = handles.D.edf;
%     %% calculate predictive matrix
    disp('Generating the Matrix of Predictors...');
    predictorMatrix = generateSWDPredictors_mj( EDF );
    
    %% load SVM classifier
    % savepath = strcat(cookieMonster, '/jonesLab_data/sleep_and_seizures/EEG_data/RQ/SWDClassificationData/');
    savepath = '/Users/Shared/JonesLab_Git/joneslab_code/Projects/sleep_and_seizures/EEG_Analysis/seizureIdentification/SVMClassifier_SWDs/detectSWDs GUI/'
    load(strcat(savepath, 'polynomialSWDSVMClassifier.mat'));

    %% classify the SWDs
    [yfit, ~] = trainedClassifier.predictFcn(predictorMatrix); 
% %     tmpfig = figure;
% %     plot(yfit)    
% %     pause
% %     close(tmpfig)
    
    %% Package and pass to epochpicker GUI
    handles.D.SVM.predictorMatrix = predictorMatrix;
    handles.D.SVM.yfit = yfit;
%     handles.D.SVM.predictorMatrix = rand(3);
%     dum1 = rand(1,handles.D.edf.nrecs);
%     dum2 = find(dum1>0.9999);
%     dum3 = 0.*dum1;
%     dum3(dum2) = 1;
%     handles.D.SVM.yfit = dum3;
    
% %% seizure viewer
% 
% %key = 1;
% loopStart = 1;
% i = 1;
% signal = EDF.D.edf.signalMat(:, 2);                 % WHY Only THIS channel?
% sampleRate = 1024; % per 4 second epoch
% recordCmp.edfLocs = find(yfit) * sampleRate;
% figure('units', 'inch', 'pos', [1 5 20 5]);
% 
%   
% while i < length(recordCmp.edfLocs)
%     for i = loopStart:length(recordCmp.edfLocs)
%         tagStart = recordCmp.edfLocs(i);
%             if tagStart < 3 * sampleRate
%                 recordCmp.exactSeizureStart(i) = 0; % place start locations in object
%                 recordCmp.exactSeizureStop(i) = 0;
%             else
%                 beforeWin = sampleRate * 2 - 1;
%                 afterWin = sampleRate * 2 - 1;
%                 seq = (tagStart - beforeWin):(tagStart + afterWin);
%                 plot(signal(seq));
%                 minimum = min(signal(seq)) + (min(signal(seq)) * 0.1); % minimum value for the EEG seq we are looking at plus 10 percent
%                 maximum = max(signal(seq)) + (max(signal(seq)) * 0.1); % same but for max
%                 bottom = -1 * max(abs([minimum, maximum])); % calculates where the bottom of the rectange that highlights the active epoch should go
%                 top = 2 * abs(bottom); % determines the top of the rectange
%                 rectangle('Position', [sampleRate, bottom, sampleRate, top], 'EdgeColor', 'r'); % places a rectangle on the plot that indicates the orignial tag location
%                 set(gca, 'xlim', [0, sampleRate*3]);
%                 
%                 % plot some labels and information
%                 text(1100, abs(bottom) - 0.01, ['record ', num2str(i), ' of ', num2str(length(recordCmp.edfLocs))]);
%                 %text(1100, bottom + .01, ['Distance from separator: ',
%                 %num2str(abs(distanceFromClassifierThresh(i, 1)))]); % not
%                 %sure why this is not represetitive of the distance from
%                 %the threshold, but its not..
% 
%                 %waitforbuttonpress;
%                 [a, ~, key] = ginput(2); % code that allows user to input 2 clicks for start and stop
%                 recordCmp.exactSeizureStart(i) = round(a(1,1)) + (tagStart - beforeWin); % place start locations in object
%                 recordCmp.exactSeizureStop(i) = round(a(2,1)) + (tagStart - beforeWin); % place stop locations in object
% 
% 
%                 % go backwards if you press the spacebar 2x
%                 if key ~= 1
%                     loopStart = i - 1;
%                     break
%                 end
%             end
%     end
% end
% close all
    
    
    % Open epochpicker GUI
        guidata(handles.EEG_GUI, handles);
        assignin('base', 'handles', handles);
%         
%         disp('Saving TEMP data after feature computation and SVM classification')
%         save([savepath 'DELETE_ME_TempOutput.mat'], 'handles', '-v7.3')
        
%% Start epochpicker GUI
            disp('Initializing epochpicker GUI')
            epochpicker_output = epochpicker(handles.EEG_GUI, [], handles);
%% --
        handles.D = epochpicker_output.D;
        disp('Done with epochpicker')
            
        guidata(handles.EEG_GUI, handles);
        assignin('base', 'handles', handles);
    
    


% --------------------------------------------------------------------
function Load_Training_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Training_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Train_SVM_Callback(hObject, eventdata, handles)
% hObject    handle to Train_SVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_Trained_SVM_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Trained_SVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_Trained_SVM_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Trained_SVM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
