function varargout = mesoRD_unitconversionGUI(varargin)
% MESORD_UNITCONVERSIONGUI MATLAB code for mesoRD_unitconversionGUI.fig
%      MESORD_UNITCONVERSIONGUI, by itself, creates a new MESORD_UNITCONVERSIONGUI or raises the existing
%      singleton*.
%
%      H = MESORD_UNITCONVERSIONGUI returns the handle to a new MESORD_UNITCONVERSIONGUI or the handle to
%      the existing singleton*.
%
%      MESORD_UNITCONVERSIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MESORD_UNITCONVERSIONGUI.M with the given input arguments.
%
%      MESORD_UNITCONVERSIONGUI('Property','Value',...) creates a new MESORD_UNITCONVERSIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mesoRD_unitconversionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mesoRD_unitconversionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mesoRD_unitconversionGUI

% Last Modified by GUIDE v2.5 12-Dec-2014 15:59:12

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesoRD_unitconversionGUI.m, sub-GUI in the SMeagol package
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
                   'gui_OpeningFcn', @mesoRD_unitconversionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mesoRD_unitconversionGUI_OutputFcn, ...
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
% --- Executes just before mesoRD_unitconversionGUI is made visible.
function mesoRD_unitconversionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mesoRD_unitconversionGUI (see VARARGIN)

% Choose default command line output for mesoRD_unitconversionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% get the text for the mesoRD input file
trj=varargin{1};                % input trj object
trj=rmfield(trj,{'speciesNames','degradedName'});
trjStrings=SM_opt2str(trj);  % convert struct to string cell vector
fileDisplay=sprintf('parsed parameters:\n--------------------');
for k=1:length(trjStrings)
    fileDisplay=sprintf('%s\n%s',fileDisplay,trjStrings{k});    
end
set(handles.text_MesoRD_log,'String',fileDisplay);
set(handles.figure1,'Name','MesoRD unit conversion');

% UIWAIT makes mesoRD_unitconversionGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);
% --- Outputs from this function are returned to the command line.
function varargout = mesoRD_unitconversionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;

varargout{1}=str2double(get(handles.set_length_factor,'String'));
varargout{2}=str2double(get(handles.set_time_factor,'String'));

% The figure can be deleted now
delete(handles.figure1);
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end
function set_length_factor_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function set_length_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_length_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','1e7'); % default value: converts cm to nm.

function set_time_factor_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function set_time_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_time_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on button press in close_GUI.
function close_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to close_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(handles.figure1);
else
    % The GUI is no longer waiting, just close it
    delete(handles.figure1);
end
