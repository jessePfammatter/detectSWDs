function varargout = epochpicker(varargin)
% EPOCHPICKER MATLAB code for epochpicker.fig
%      EPOCHPICKER, by itself, creates a new EPOCHPICKER or raises the existing
%      singleton*.
%
%      H = EPOCHPICKER returns the handle to a new EPOCHPICKER or the handle to
%      the existing singleton*.
%
%      EPOCHPICKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EPOCHPICKER.M with the given input arguments.
%
%      EPOCHPICKER('Property','Value',...) creates a new EPOCHPICKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before epochpicker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to epochpicker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help epochpicker

% Last Modified by GUIDE v2.5 13-Oct-2016 21:08:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @epochpicker_OpeningFcn, ...
                   'gui_OutputFcn',  @epochpicker_OutputFcn, ...
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


% --- Executes just before epochpicker is made visible.
function epochpicker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to epochpicker (see VARARGIN)

% Choose default command line output for epochpicker
handles.output = hObject;

% Extract the passed data and package into a LOCAL handles struct
dummy = varargin{3};
handles.D = dummy.D;
handles.D.status = 'epochpicker opened';
% signal = handles.D.edf.signalMat(:, 2);     	% WHY Only THIS channel?
signal = handles.D.edf.signalMat;             	% WHY Only THIS channel?
norm_signal = -0.1 .* signal ./ repmat(std(signal), size(signal,1), 1);
norm_signal = norm_signal - repmat([1:size(norm_signal, 2)], size(norm_signal,1), 1);
norm_signal = norm_signal - mean(norm_signal(:));
signal = norm_signal;
epochdur = handles.D.edf.epochdur;
fs = handles.D.edf.fs;
npts = 1024 .* epochdur;    %%% SOMETHING IS WRONG HERE, fs, epochdur, recdur, data_record_duration and npts DO NOT ADD UP PROPERLY in handles.D.edf.
                            %%% This should NOT have to be hardcoded like
                            %%% this. It is dangerous and will lead to
                            %%% errors later.
nepochs = handles.D.edf.nrecs./epochdur-1;
epochnums = [1:nepochs];
recordCmp.edfLocs = find(handles.D.SVM.yfit) * npts;
hit_epochs = find(handles.D.SVM.yfit);
curr_hit = 0;

    % beforeWin = sampleRate * 4 - 1;
% afterWin = sampleRate * 4;
preepochs   = 1;
postepochs  = 2;
winpts      = [-preepochs.*npts : postepochs.*npts-1];

% Draw the initial display
axes(handles.Epoch_Axes)
    curr_epoch = 1 + preepochs; % first displayable epoch
	
    epochstartloc = (curr_epoch-1).*npts+1;

        seq = epochstartloc + winpts;
        timeseg = seq./fs./epochdur;
        plot(timeseg, signal(seq, :));
        mi = min(signal(seq,:)) + (min(signal(seq,:)) * 0.1); % minimum value for the EEG seq we are looking at plus 10 percent
        mx = max(signal(seq,:)) + (max(signal(seq,:)) * 0.1); % same but for max
        bottom = -1 * max(abs([mi, mx])); % calculates where the bottom of the rectange that highlights the active epoch should go
        top = 2 * abs(bottom); % determines the top of the rectange
        rectangle('Position', [epochstartloc./npts.*epochdur, bottom, epochdur, top], 'EdgeColor', 'r'); % places a rectangle on the plot that indicates the orignial tag location
        set(gca, 'xlim', timeseg([1, end]), 'ylim', handles.D.edf.nchans.*[-0.6 0.6]);

        % plot some labels and information
        text(epochstartloc./fs./epochdur, 0.9.*abs(bottom), ['Epoch #', num2str(curr_epoch), ' of ', num2str(nepochs)], 'fontsize', 18);

    xlabel('Time, sec')
    ylabel('Normalized Amplitude')    
    hold off

axes(handles.Hit_Axes)
    hit_status = ismember(curr_epoch, hit_epochs);
    plot(epochnums(:), handles.D.SVM.yfit(:), 'k-'); hold on
    plot(curr_epoch, hit_status, 'ro', 'markerfacecolor', 'r')
    set(gca, 'xlim', epochnums([1,end]), 'ylim', [0 1.1])
    box off
    xlabel('Epoch #')
    ylabel('Hits')
    hold off

    % Package into guidata
    handles.D.SVM.temp.signal               = signal;
    handles.D.SVM.temp.npts                 = npts;
    handles.D.SVM.temp.hit_epochs           = hit_epochs;    
    handles.D.SVM.temp.curr_hit             = curr_hit;
    handles.D.SVM.temp.curr_epoch         	= curr_epoch;
    handles.D.SVM.temp.recordCmp.edfLocs    = recordCmp.edfLocs;
    handles.D.SVM.temp.preepochs            = preepochs;
    handles.D.SVM.temp.postepochs           = postepochs;
    handles.D.SVM.temp.winpts               = winpts;
    handles.D.SVM.temp.nepochs              = nepochs;
    handles.D.SVM.temp.epochnums            = epochnums;
    handles.D.SVM.temp.epochdur             = epochdur;
    handles.D.SVM.temp.fs                   = fs;


% Update handles structure
guidata(hObject, handles);
chirpsound_mj

% UIWAIT makes epochpicker wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = epochpicker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
handles.D.status = 'epochpicker finished';
varargout{1} = handles;




% --- Executes on slider movement.
function Epoch_Slider_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
    nepochs = handles.D.edf.nsecs./handles.D.edf.epochdur;
    curr_epoch = round(get(hObject,'Value').*nepochs);
        if curr_epoch < 1; 
            curr_epoch = 1;
        elseif curr_epoch > nepochs 
            curr_epoch = nepochs-1;
        end
    handles.D.SVM.temp.curr_epoch = curr_epoch;    
% Update handles structure and redraw
guidata(hObject, handles);
Update_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function Epoch_Slider_CreateFcn(hObject, eventdata, handles)
%     nepochs = handles.D.edf.nsecs./handles.D.edf.epochdur;
%     mi = 1./nepochs;
%     mx = 100./nepochs;
%     set(hObject,'SliderStep',[mi mx]);




function Update_Callback(hObject, eventdata, handles)
    signal              = handles.D.SVM.temp.signal;
    npts                = handles.D.SVM.temp.npts;
    curr_hot            = handles.D.SVM.temp.curr_hit;
    curr_epoch          = handles.D.SVM.temp.curr_epoch;
    hit_epochs          = handles.D.SVM.temp.hit_epochs;
    recordCmp.edfLocs   = handles.D.SVM.temp.recordCmp.edfLocs;
    preepochs           = handles.D.SVM.temp.preepochs;
    postepochs          = handles.D.SVM.temp.postepochs;
    winpts              = handles.D.SVM.temp.winpts;
    nepochs             = handles.D.SVM.temp.nepochs;
    epochnums           = handles.D.SVM.temp.epochnums;
    epochdur           = handles.D.SVM.temp.epochdur;
    fs                 = handles.D.SVM.temp.fs;

     	if curr_epoch - preepochs < 1; % Prevent from running off the edges
            curr_epoch = curr_epoch + preepochs;
        elseif curr_epoch + postepochs > nepochs 
            curr_epoch = nepochs-postepochs;
        end
    epochstartloc = (curr_epoch-1).*npts+1;

    seq = epochstartloc + winpts;
    timeseg = seq./fs./epochdur;
    
    axes(handles.Epoch_Axes)
    plot(timeseg, signal(seq, :));
    mi = min(signal(seq,:)) + (min(signal(seq,:)) * 0.1); % minimum value for the EEG seq we are looking at plus 10 percent
    mx = max(signal(seq,:)) + (max(signal(seq,:)) * 0.1); % same but for max
    bottom = -1 * max(abs([mi, mx])); % calculates where the bottom of the rectange that highlights the active epoch should go
    top = 2 * abs(bottom); % determines the top of the rectange
    rectangle('Position', [epochstartloc./npts.*epochdur, bottom, epochdur, top], 'EdgeColor', 'r'); % places a rectangle on the plot that indicates the orignial tag location
    set(gca, 'xlim', timeseg([1, end]), 'ylim', handles.D.edf.nchans.*[-0.6 0.6]);

    % plot some labels and information
    text(epochstartloc./fs./epochdur, 0.9.*abs(bottom), ['Epoch #', num2str(curr_epoch), ' of ', num2str(nepochs)], 'fontsize', 18);

    xlabel('Time, sec')
    ylabel('Normalized Amplitude')  
    hold off

axes(handles.Hit_Axes)
    hit_status = ismember(curr_epoch, hit_epochs);
    plot(epochnums(:), handles.D.SVM.yfit(:), 'k-'); hold on
    plot(curr_epoch, hit_status, 'ro', 'markerfacecolor', 'r')
    set(gca, 'xlim', epochnums([1,end]), 'ylim', [0 1.1])
    box off
    xlabel('Epoch #')
    ylabel('Hits')
    hold off
    
    
    
    



function Epoch_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to Epoch_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Epoch_Edit as text
%        str2double(get(hObject,'String')) returns contents of Epoch_Edit as a double


% --- Executes during object creation, after setting all properties.
function Epoch_Edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Epoch_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Epoch_Jump_Back_Button.
function Epoch_Jump_Back_Button_Callback(hObject, eventdata, handles)
    handles.D.SVM.temp.curr_epoch  = handles.D.SVM.temp.curr_epoch - 1;
        if handles.D.SVM.temp.curr_epoch < 2; 
            handles.D.SVM.temp.curr_epoch = 2; 
        elseif handles.D.SVM.temp.curr_epoch > handles.D.SVM.temp.nepochs; 
            handles.D.SVM.temp.curr_epoch = handles.D.SVM.temp.nepochs; 
        end
    set(handles.Epoch_Edit, 'str', num2str(handles.D.SVM.temp.curr_epoch) )
    set(handles.Epoch_Slider, 'val', round(handles.D.SVM.temp.curr_epoch./handles.D.SVM.temp.nepochs) ) 
    % Update handles structure and redraw
    guidata(hObject, handles);
    Update_Callback(hObject, eventdata, handles)


% --- Executes on button press in Epoch_Jump_Forward_Button.
function Epoch_Jump_Forward_Button_Callback(hObject, eventdata, handles)
    handles.D.SVM.temp.curr_epoch  = handles.D.SVM.temp.curr_epoch + 1;
        if handles.D.SVM.temp.curr_epoch < 2; 
            handles.D.SVM.temp.curr_epoch = 2; 
        elseif handles.D.SVM.temp.curr_epoch > handles.D.SVM.temp.nepochs; 
            handles.D.SVM.temp.curr_epoch = handles.D.SVM.temp.nepochs; 
        end
    set(handles.Epoch_Edit, 'str', num2str(handles.D.SVM.temp.curr_epoch) )
    set(handles.Epoch_Slider, 'val', handles.D.SVM.temp.curr_epoch./handles.D.SVM.temp.nepochs ) 
    % Update handles structure and redraw
    guidata(hObject, handles);
    Update_Callback(hObject, eventdata, handles)


% --- Executes on button press in Hit_Jump_Back_Button.
function Hit_Jump_Back_Button_Callback(hObject, eventdata, handles)
    handles.D.SVM.temp.curr_hit     = handles.D.SVM.temp.curr_hit - 1;
        if handles.D.SVM.temp.curr_hit < 1; 
            handles.D.SVM.temp.curr_hit = 1; 
        elseif handles.D.SVM.temp.curr_hit > length(handles.D.SVM.temp.hit_epochs); 
            handles.D.SVM.temp.curr_hit = length(handles.D.SVM.temp.hit_epochs); 
        end
    handles.D.SVM.temp.curr_epoch   = handles.D.SVM.temp.hit_epochs(handles.D.SVM.temp.curr_hit);
    set(handles.Epoch_Edit, 'str', num2str(handles.D.SVM.temp.curr_epoch) )
    set(handles.Epoch_Slider, 'val', round(handles.D.SVM.temp.curr_epoch./handles.D.SVM.temp.nepochs) ) 
    % Update handles structure and redraw
    guidata(hObject, handles);
    Update_Callback(hObject, eventdata, handles)


% --- Executes on button press in Hit_Jump_Forward_Button.
function Hit_Jump_Forward_Button_Callback(hObject, eventdata, handles)
    handles.D.SVM.temp.curr_hit     = handles.D.SVM.temp.curr_hit + 1;
        if handles.D.SVM.temp.curr_hit < 1; 
            handles.D.SVM.temp.curr_hit = 1; 
        elseif handles.D.SVM.temp.curr_hit > length(handles.D.SVM.temp.hit_epochs); 
            handles.D.SVM.temp.curr_hit = length(handles.D.SVM.temp.hit_epochs); 
        end
    handles.D.SVM.temp.curr_epoch   = handles.D.SVM.temp.hit_epochs(handles.D.SVM.temp.curr_hit);
    set(handles.Epoch_Edit, 'str', num2str(handles.D.SVM.temp.curr_epoch) )
    set(handles.Epoch_Slider, 'val', round(handles.D.SVM.temp.curr_epoch./handles.D.SVM.temp.nepochs) ) 
    % Update handles structure and redraw
    guidata(hObject, handles);
    Update_Callback(hObject, eventdata, handles)


% --- Executes on button press in Mark_Hits_No_Text.
function Mark_Hits_No_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Mark_Hits_No_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Mark_Hits_Yes_Button.
function Mark_Hits_Yes_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Mark_Hits_Yes_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
