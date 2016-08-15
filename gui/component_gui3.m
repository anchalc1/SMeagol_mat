function varargout = component_gui3(varargin)
% opt=component_gui3(componentName,inputOpt)
% 
% opens a matlab GUI to display and edit an options structure associated
% with a Palantir component.
% input:
% componentName     : name of the palantir component
% inputOpt          : initial opt struct for componentName. If not given,
% component_gui3 starts with the default options given by function Name.
%
% ML 2014-12-17

% COMPONENT_GUI3 M-file for component_gui3.fig
%      COMPONENT_GUI3, by itself, creates a new COMPONENT_GUI3 or raises the existing
%      singleton*.
%
%      H = COMPONENT_GUI3 returns the handle to a new COMPONENT_GUI3 or the handle to
%      the existing singleton*.
%
%      COMPONENT_GUI3('Property','Value',...) creates a new COMPONENT_GUI3 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to component_gui3_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      COMPONENT_GUI3('CALLBACK') and COMPONENT_GUI3('CALLBACK',hObject,...) call the
%      local function named CALLBACK in COMPONENT_GUI3.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help component_gui3

% Last Modified by GUIDE v2.5 18-Dec-2014 14:02:38

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% component_gui3.m, sub-GUI in the SMeagol package
% ========================================================================= 
% Copyright (C) 2015 Martin Lindén and Johan Elf
% 
% E-mail: bmelinden@gmail.com, johan.elf@gmail.com
% =========================================================================
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or any
% later version.  This program is distributed in the hope that it will
% be useful, but WITHOUT ANY WARRANTY; without even the implied
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
% the GNU General Public License for more details.
% 
% Additional permission under GNU GPL version 3 section 7
%  
% If you modify this Program, or any covered work, by linking or
% combining it with Matlab or any Matlab toolbox, the licensors of this
% Program grant you additional permission to convey the resulting work.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @component_gui3_OpeningFcn, ...
                   'gui_OutputFcn',  @component_gui3_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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
% --- Executes just before component_gui3 is made visible.
function component_gui3_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to component_gui1 (see VARARGIN)

% Choose default command line output for component_gui1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes component_gui1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% get GUI data structure
h=guidata(hObject);

% get function name and put it in the title bar
if(nargin<=3) % no extra input, manually select a function and use default options
    [SMfile,SMpath]=uigetfile(...
    {'SM_*.m','Palantir function (SM_*.m)';
    '*.m','matlab functiojn (*.m)'},'Select Palantir file for options editing');
    [~,functionName,~]=fileparts(fullfile(SMpath,SMfile));
    inputOpt=eval([functionName '()']);
elseif(nargin==4) % function name given, but no input options
    functionName=varargin{1};
    inputOpt=eval([functionName '()']); % use default as input
elseif(nargin==5) % both function name and input parameters given
    functionName=varargin{1};
    inputOpt=varargin{2};
elseif(nargin>5)   % weird
    error('component_gui1 called with too many input argument.')
end
set(hObject,'Name',['Options fields for ' functionName])

defaultOpt=eval([functionName '()']); % get default options from the function itself

% create bare-bone options structure for output

% check input options compatible with default options
currentOpt=struct;
defaultFields=fieldnames(defaultOpt);
for k=1:length(defaultFields)
    if(~isfield(inputOpt,defaultFields(k)))
        error('component_gui1 called with incompatible function name and options structure')    
    else % transfer the field to current options opbject
        currentOpt.(defaultFields{k})=inputOpt.(defaultFields{k});
    end
end

% make the edit box be left-aligned and have multiple lines (Max-Min>1)
set(h.optionTextEdit,'Min',0);
set(h.optionTextEdit,'Max',2);
set(h.optionTextEdit,'HorizontalAlignment','left');

% add fields to guidata
h.functionName=functionName;
h.helptext=help(functionName);
h.inputOpt  =  inputOpt;
h.defaultOpt=defaultOpt;
h.currentOpt=currentOpt;

% insert text representation of current options
set(h.optionTextEdit,'String',SM_opt2str(h.currentOpt));

% save modified GUI data
guidata(hObject,h);
% --- Outputs from this function are returned to the command line.
function varargout = component_gui3_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiwait(handles.figure1); % wait here until closing time
h=guidata(hObject);
varargout{1}=h.currentOpt;
delete(h.figure1);

% --- Executes on button press in restoreDefault.
function restoreDefault_Callback(hObject, eventdata, handles)
% hObject    handle to restoreDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=guidata(hObject);
h.currentOpt=h.defaultOpt;                          % new current options
set(h.optionTextEdit,'String',SM_opt2str(h.currentOpt));% update text display
guidata(hObject,h);                                 % store GUI data
% --- Executes on button press in saveAndClose.
function saveAndClose_Callback(hObject, eventdata, handles)
% hObject    handle to saveAndClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1); % resume in component_gui3_OutputFcn
% --- Executes on button press in showHelp.
function showHelp_Callback(hObject, eventdata, handles)
% hObject    handle to showHelp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get GUI data structure
h=guidata(hObject);
helpdlg(help(h.functionName),h.functionName);
% --- Executes on button press in restoreInput.
function restoreInput_Callback(hObject, eventdata, handles)
% hObject    handle to restoreInput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=guidata(hObject);
h.currentOpt=h.inputOpt; % change current options
set(h.optionTextEdit,'String',SM_opt2str(h.currentOpt));% update text display
guidata(hObject,h); % store new settings
function optionTextEdit_Callback(hObject, eventdata, handles)
% hObject    handle to optionTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of optionTextEdit as text
%        str2double(get(hObject,'String')) returns contents of optionTextEdit as a double

set(handles.saveAndClose,'Enable','off');
set(handles.saveAndClose,'String','Evaluate before closing.');
guidata(hObject,handles); % store changes
% --- Executes during object creation, after setting all properties.

function optionTextEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to optionTextEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in addFile.
function addFile_Callback(hObject, eventdata, handles)
% hObject    handle to addFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[newFile,newPath]=uigetfile(...
       {'*.m;*.fig;*.mat', 'MATLAB File (*.m, *.fig, *.mat)'; ...
        '*.tif;*tiff',  'tif file (*.tif, *.tiff)'; ...
        '*.m',  'MATLAB Code (*.m)'; ...
        '*.mat','MAT-files (*.mat)'; ...
        '*',  'All Files (*)'},'Select file to add as a parameter');
% if file selection is cancelled newFile and newPath =0 (not char)
if(ischar(newFile) && ischar(newPath))
    parString=get(handles.optionTextEdit,'String');
    
    % select variable to put the new path in
    oldpar=fieldnames(handles.currentOpt);
    [selection,OK]=listdlg('ListString',{oldpar{:},'new_parameter'},...
                           'SelectionMode','single',...
                           'InitialValue',length(oldpar)+1,...
                           'name','parameter name');
    
    if(~OK)
        return; % something went wrong
    elseif(OK && selection > length(oldpar))
        % then add a new parameter at the last line
        parString{end+1}=['newfile' int2str(length(parString)) '=' char(39) fullfile(newPath,newFile) char(39) ';'];
        set(handles.optionTextEdit,'String',parString);
    else % put the selected file in an existing parameter
        parString{end+1}=[oldpar{selection} '=' char(39) fullfile(newPath,newFile) char(39) ';'];
        set(handles.optionTextEdit,'String',parString);
    end
    % force reevaluation
    optionTextEdit_Callback(handles.saveAndClose,[],handles)
end
% --- Executes on button press in evaluateOptions.
function evaluateOptions_Callback(hObject, eventdata, handles)
% hObject    handle to evaluateOptions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newString=get(handles.optionTextEdit,'String');
[newOpt,isOK]=SM_str2opt(newString);

if(~isOK)
   errordlg('the text could not be converted. Reverting to latest correct text')
else
    % check options consistency with default object?
    handles.currentOpt=newOpt;
    
    % maybe keep text format as is as long as it is readable?
    % handles.optionTextEdit.String=SM_opt2str(handles.currentOpt);
end
set(handles.optionTextEdit,'String',SM_opt2str(handles.currentOpt));
set(handles.saveAndClose,'Enable','on');
set(handles.saveAndClose,'String','close');
guidata(hObject,handles); % store changes
