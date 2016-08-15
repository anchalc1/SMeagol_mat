function varargout = SM_mesoRD_gui(varargin)
%SM_MESORD_GUI M-file for SM_mesoRD_gui.fig
%      SM_MESORD_GUI, by itself, creates a new SM_MESORD_GUI or raises the existing
%      singleton*.
%
%      H = SM_MESORD_GUI returns the handle to a new SM_MESORD_GUI or the handle to
%      the existing singleton*.
%
%      SM_MESORD_GUI('Property','Value',...) creates a new SM_MESORD_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SM_mesoRD_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SM_MESORD_GUI('CALLBACK') and SM_MESORD_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SM_MESORD_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SM_mesoRD_gui

% Last Modified by GUIDE v2.5 10-Sep-2015 15:09:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SM_mesoRD_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @SM_mesoRD_gui_OutputFcn, ...
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


% --- Executes just before SM_mesoRD_gui is made visible.
function SM_mesoRD_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SM_mesoRD_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SM_mesoRD_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% Choose default command line output for SMeagol_gui

%% initialize options structure and reset GUI values
h=guidata(hObject);

% remember path the file was called from
h.gui_start_path=pwd;
guidata(hObject,h);

% reset fields to default values
restoreDefaultParameters(hObject);
% --- Outputs from this function are returned to the command line.
function varargout = SM_mesoRD_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ff=msgbox('Detecting MesoRD with tracking features.')
R=SM_mesoRD_tracking_test();
try
    delete(ff);
catch
end
if(R)
    uiwait(handles.figure1); % wait here until closing time
else
    rr=errordlg('MesoRD tracking test failed. See SM_mesoRD_tracking_test and SM_find_mesoRD_path.');
end
% try to extract location of mesoRD logfile
h=guidata(hObject);
logfile=(fullfile(get(h.target_folder_txt,'String'),'molecule_tracking_log.txt'));
if(exist(logfile,'file'))
    h.output=logfile;
else
    h.output='';
end
varargout{1}=h.output;
delete(h.figure1);

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%% internal functions
function restoreDefaultParameters(hObject)
% restore all fields to default values, sort of a quick GUI restart
h=guidata(hObject);
set(h.SBML_model_txt,'String','-');
set(h.target_folder_txt,'String',h.gui_start_path);
set(h.voxel_size_edit,'String','0.010 um');
set(h.tracking_timestep_edit,'String','0.005');
set(h.simulation_time_edit,'String','5');
set(h.tracking_init_time_edit,'String','0');

set(h.misc_flags_txt,'String',' -i 1 -I 50 -c -1 -C -1 -E -g -K');
%set(h.misc_flags_txt,'String',' -i 1 -I 50 -c -1 -C -1 -E -p -g -K');
% -p does not work on windows

guidata(hObject,h);
%% callback and creation functions
% --- Executes on button press in SBML_model_button.
function SBML_model_button_Callback(hObject, eventdata, handles)
% hObject    handle to SBML_model_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% select start path: 1) existing SBML file 2) target folder 3) gui start
% folder:
h=guidata(hObject);
[SBMLpath,~,~]=fileparts(get(h.SBML_model_txt,'String'));
[targetPath,~,~]=fileparts(fullfile(get(h.target_folder_txt,'String'),'.'));

if(    ~isempty(SBMLpath) && exist(SBMLpath,'dir'))
    cd(SBMLpath)
elseif(~isempty(targetPath) && exist(targetPath,'dir'))
    cd(targetPath)
else
    cd(h.gui_start_path)
end
[modelfile,modelpath]=uigetfile(...
    {'*.xml','MesoRD SBML file';
    '*.txt','text file (.txt)';'*','All files'},'select SBML model file');
cd(h.gui_start_path);
set(h.SBML_model_txt,'String',fullfile(modelpath,modelfile));
guidata(hObject,h);
% --- Executes during object creation, after setting all properties.
function SBML_model_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SBML_model_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in target_folder_button.
function target_folder_button_Callback(hObject, eventdata, handles)
% hObject    handle to target_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h=guidata(hObject);
[SBMLpath,~,~]=fileparts(get(h.SBML_model_txt,'String'));
[targetPath,~,~]=fileparts(fullfile(get(h.target_folder_txt,'String'),'.'));
% select start path: 1) existing target folder 2) SBML file folder 3) gui
% start folder:
if(~isempty(targetPath) && exist(targetPath,'dir'))
    cd(targetPath)
elseif(    ~isempty(SBMLpath) && exist(SBMLpath,'dir'))
    cd(SBMLpath)
else
    cd(h.gui_start_path)
end
targetdir=uigetdir('','select simulation target folder');
cd(h.gui_start_path);
set(h.target_folder_txt,'String',targetdir);
guidata(hObject,h);
% --- Executes during object creation, after setting all properties.
function target_folder_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_folder_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function text100_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function text5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in run_simulation_button.
function run_simulation_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_simulation_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


h=guidata(hObject);
% go to target path
cd(get(h.target_folder_txt,'String'));
if(~exist(get(h.SBML_model_txt,'String'),'file'))
    errordlg('Model file not specified or does not exist');
    return;
end

mesoRDoptions=[ '"' get(h.SBML_model_txt,'String') '" ' ... % hack to deal with pathnames with spaces
                get(h.misc_flags_txt,'String') ' ' ...
                ' -q ' get(h.voxel_size_edit,'String')...
                ' -w ' get(h.tracking_timestep_edit,'String')...
                ' -t ' get(h.simulation_time_edit,'String')...
                ' -x ' get(h.tracking_init_time_edit,'String')];%' & '];

b=questdlg('MesoRD simulations can take a long time to complete. Continue?','Yes','No');
if(~strcmp(b,'Yes'))
    return;
end
hd=helpdlg('MesoRD output in the commandline','MesoRD options');

sysflag=SM_runMesoRD(mesoRDoptions);
try
    delete(hd);
catch
end

cd(h.gui_start_path)
% --- Executes during object creation, after setting all properties.
function run_simulation_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_simulation_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not crgeated until after all CreateFcns called
% --- Executes on button press in close_gui_button.
function close_gui_button_Callback(hObject, eventdata, handles)
% hObject    handle to close_gui_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%warning('GUI closure and output not implemented')
uiresume(handles.figure1);
%delete(handles.figure1);
% --- Executes during object creation, after setting all properties.
function close_gui_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to close_gui_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function SBML_model_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SBML_model_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function target_folder_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to target_folder_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
function voxel_size_edit_Callback(hObject, eventdata, handles)
% hObject    handle to voxel_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of voxel_size_edit as text
%        str2double(get(hObject,'String')) returns contents of voxel_size_edit as a double
% --- Executes during object creation, after setting all properties.
function voxel_size_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to voxel_size_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tracking_timestep_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tracking_timestep_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tracking_timestep_edit as text
%        str2double(get(hObject,'String')) returns contents of tracking_timestep_edit as a double
% --- Executes during object creation, after setting all properties.
function tracking_timestep_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracking_timestep_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function tracking_init_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tracking_init_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of tracking_init_time_edit as text
%        str2double(get(hObject,'String')) returns contents of tracking_init_time_edit as a double
% --- Executes during object creation, after setting all properties.
function tracking_init_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracking_init_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function misc_flags_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to misc_flags_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in mesord_help_button.
function mesord_help_button_Callback(hObject, eventdata, handles)
% hObject    handle to mesord_help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hh=helpdlg('MesoRD help text is shown in the commandline','MesoRD options');
system('mesord -h')
pause(3);
try
    delete(hh)
catch
end
% --- Executes during object creation, after setting all properties.
function mesord_help_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mesord_help_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes on button press in restore_defaults_button.
function restore_defaults_button_Callback(hObject, eventdata, handles)
% hObject    handle to restore_defaults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
restoreDefaultParameters(hObject);
% --- Executes during object creation, after setting all properties.
function restore_defaults_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to restore_defaults_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% --- Executes during object creation, after setting all properties.
function text9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
function simulation_time_edit_Callback(hObject, eventdata, handles)
% hObject    handle to simulation_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of simulation_time_edit as text
%        str2double(get(hObject,'String')) returns contents of simulation_time_edit as a double
% --- Executes during object creation, after setting all properties.
function simulation_time_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simulation_time_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%uiresume(handles.figure1);
delete(hObject);
% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
try
    delete(handles.figure1);
catch me
end
