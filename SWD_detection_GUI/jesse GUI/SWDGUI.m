function varargout = SWDGUI(varargin)
% SWDGUI MATLAB code for SWDGUI.fig
% SWDGUI allows a user to inspect SWD-like EEG events and decide whether or
% not they are true SWD events. This GUI works by loading a struct 'D' that
% contains data, and zz. 
%
% D.data is an data structure that contains all the information to run the
% GUI. An example data set is provided at
% '~/Dropbox/research/jonesLab/jonesLab_git/jonesLab_code/Projects/sleep_and_seizures/EEG_Analysis/seizureIdentification/SWD_detection/SWD_detection_GUI/SWDtest/SWDGUI_exampleData.mat'
%
% D.zz is the starting point index. D.zz updates as user works through data
% set and allows user to pick up where they left off upon saving, closing,
% and repoinging the GUI.
%
% To run this file, open the SWDGUI.m script, and choose an appropriate
% data set. Complete scoreing of SWDs (and save progress along the way).
% Save progress and close. Data from scoring is saved in D.data.responses
% and D.data.comments. Data in 'SWDGUI_exampleData.mat' were generated from
% the SVM marked epochs where SWD start and stops were identified with
% 'batch_detectSWDs.m'. Start and stop times were proofed by JP with
% 'proofSWD.m' and the data set was compiled on 4/24/17. 
%
% JP 2017

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SWDGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SWDGUI_OutputFcn, ...
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

% --- Executes just before SWDGUI is made visible.
function SWDGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SWDGUI (see VARARGIN)

% load data

[handles.filename, handles.pathname] = uigetfile();
load(strcat(handles.pathname, handles.filename));

% assign handles.D
handles.D = D;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.output = hObject;
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = SWDGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in accept.
function accept_Callback(hObject, eventdata, handles)
% hObject    handle to accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% mark as accept
D.data(D.zz).responses = 2;

%plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);


% --- Executes on button press in decline.
function decline_Callback(hObject, eventdata, handles)
% hObject    handle to decline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% mark as decline
D.data(D.zz).responses = 1;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);



% --- Executes on button press in previous.
function previous_Callback(hObject, eventdata, handles)
% hObject    handle to previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% go to previous event
D.zz = D.zz - 1;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);

% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% go to next event
D.zz = D.zz + 1;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);

function whichSWD_Callback(hObject, eventdata, handles)
% hObject    handle to whichSWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% jump to event entered
D.zz = str2num(get(hObject, 'String'));

% add warning if you enter something invalid
if  D.zz > size(D.data, 2) 
   warndlg(strcat({'Warning: Entry must be a number between 1 and '}, num2str(size(D.data, 2)), '!'));
   pause
end

if  D.zz < 1
   warndlg(strcat({'Warning: Entry must be a number between 1 and '}, num2str(size(D.data, 2)), '!'));
   pause
end

if  isempty(D.zz)
   warndlg(strcat({'Warning: Entry must be a number between 1 and '}, num2str(size(D.data, 2)), '!'));
   pause
end


% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichSWD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichSWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function comments_Callback(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% assign handles.D
D = handles.D;

% store comments
D.data(D.zz).comments = get(hObject, 'String');

% update GUI
handles.D = D;
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function comments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to whichSWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichSWD as text
%        str2double(get(hObject,'String')) returns contents of whichSWD as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichSWD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function plottttttttt(hObject, eventdata, handles, D)
% plot the signal in the first window

% add some text here to jump to accepted or rejected SWDs by incrementing
% or decrementing D.zz

% XXX


axes(handles.eventViewingPlot);
plot(D.data(D.zz).signalClips(:, 1) + 200);
hold on;
plot(D.data(D.zz).signalClips(:, 2) - 200);
xlim([0, 2048]);
ylim([-1000, 1000]);
%text(100, 800, strcat({'Event '}, num2str(D.zz)), 'Color', 'k'); % 256 is at 1 second (remember I change the scale below)
bufferbuffer = 30;
if D.data(D.zz).responses == 0
    patch([1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer], [-1000, -1000, 1000, 1000], 'y', 'FaceColor', 'y', 'FaceAlpha', .25, 'EdgeColor', 'none');
elseif D.data(D.zz).responses == 1
    patch([1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer], [-1000, -1000, 1000, 1000], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none');
elseif D.data(D.zz).responses == 2
    patch([1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 + (.5*D.data(D.zz).seizureDuration) + bufferbuffer, 1024 - (.5*D.data(D.zz).seizureDuration) - bufferbuffer], [-1000, -1000, 1000, 1000], 'g', 'FaceColor', 'g', 'FaceAlpha', .25, 'EdgeColor', 'none');
end
hold off;
ylabel('EEG Amplitude');

% put xtick labels in seconds this is set for 256 sample rate.
% xticks([0, 256, 512, 768, 1024, 1280, 1536, 1792, 2048])
set(gca, 'xtick', [0, 256, 512, 768, 1024, 1280, 1536, 1792, 2048]) % added by MJ to avoid error
% xticklabels({'0', '1', '2', '3', '4', '5', '6', '7', '8'})
set( gca, 'xticklabels', {'0', '1', '2', '3', '4', '5', '6', '7', '8'}) % added by MJ to avoid error
xlabel('Time (s)')

% plot a numberline in second window
axes(handles.progressPlot);
for i = 1:size(D.data, 2)
    responses(i) = D.data(i).responses;
end
responseX = 1:length(responses);
plot(responseX, responses, 'k-');
hold on;
    
plot(responseX(responses == 0), responses(responses == 0), '.y', 'MarkerSize', 10);
plot(responseX(responses == 1), responses(responses == 1), '.r', 'MarkerSize', 10);
plot(responseX(responses == 2), responses(responses == 2), '.g', 'MarkerSize', 10);
plot(D.zz, 0, 'ok', 'MarkerSize', 20); % you are here marker
xlim([1 length(responses)]);
ylim([0, 2]);
xlabel('Events to Score (2 = SWD, 1 = notSWD, 0 = Need to Score)');
ylabel('Completed?');
hold off;

% set the number for display
set(handles.whichSWD, 'String', num2str(D.zz));
set(handles.text5, 'String', strcat({'of '}, num2str(size(D.data, 2))));

% display comments
if ~strcmp(D.data(D.zz).comments, 'NA')
    set(handles.comments, 'String', D.data(D.zz).comments)
else
    set(handles.comments, 'String', {' Add Comments Here: (and press Enter when finished)'})
end


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key
    case 'uparrow'
        accept_Callback(hObject, eventdata, handles)
    case 'downarrow'
        decline_Callback(hObject, eventdata, handles)
    case 'leftarrow'
        previous_Callback(hObject, eventdata, handles)
    case 'rightarrow'
        next_Callback(hObject, eventdata, handles)
    case 'u'
        unsure_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in saveProgress.
function saveProgress_Callback(hObject, eventdata, handles)
% hObject    handle to saveProgress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

D = handles.D;
save(strcat(handles.pathname, handles.filename), 'D')

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% load data

[handles.filename, handles.pathname] = uigetfile();
load(strcat(handles.pathname, handles.filename));

% assign handles.D
handles.D = D;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.output = hObject;
guidata(hObject, handles);


% --- Executes on button press in unsure.
function unsure_Callback(hObject, eventdata, handles)
% hObject    handle to unsure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% assign handles.D
D = handles.D;

% mark as decline
D.data(D.zz).responses = 0;

% plot
plottttttttt(hObject, eventdata, handles, D)

% update GUI
handles.D = D;
guidata(hObject, handles);
