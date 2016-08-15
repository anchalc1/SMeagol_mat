function varargout = SMeagol_gui(varargin)
% varargout = SMeagol_gui(varargin)
% SMEAGOL_GUI MATLAB code for SMeagol_gui.fig
%      SMEAGOL_GUI, by itself, creates a new SMEAGOL_GUI or raises the existing
%      singleton*.
%
%      H = SMEAGOL_GUI returns the handle to a new SMEAGOL_GUI or the handle to
%      the existing singleton*.
%
%      SMEAGOL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SMEAGOL_GUI.M with the given input arguments.
%
%      SMEAGOL_GUI('Property','Value',...) creates a new SMEAGOL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SMeagol_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SMeagol_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SMeagol_gui

% Last Modified by GUIDE v2.5 03-Sep-2015 09:30:46

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMeagol_gui.m, main GUI for the the SMeagol package
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
%% start of actual code

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SMeagol_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @SMeagol_gui_OutputFcn, ...
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

% --- Executes just before SMeagol_gui is made visible.
function SMeagol_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SMeagol_gui (see VARARGIN)

% Choose default command line output for SMeagol_gui
handles.output = hObject;

SM_license('SMeagol_gui');

%% initialize options structure and reset GUI values


% Palantir base path
handles.smeagol_base_path=fullfile(mfilename('fullpath'),'..','..');

% remember path the file was called from
handles.gui_start_path=pwd;

%% activation models
files=dir(fullfile(handles.smeagol_base_path,'components','photoactivation','*.m'));
models={'photoactivation model?'};
files([files.isdir]==0);
for k=1:length(files)
    [~,models{end+1}]=fileparts(files(k).name);
end    
set(handles.select_activation_model,'String',models,'value',1);
set(handles.set_activation_parameters,'Enable','off');
%% base intensity models
files=dir(fullfile(handles.smeagol_base_path,'components','baseIntensity','*.m'));
models={'base intensity model?'};
files([files.isdir]==0);
for k=1:length(files)
    [~,models{end+1}]=fileparts(files(k).name);
end    
set(handles.select_baseIntensity_model,'String',models,'value',1);
set(handles.set_baseIntensity_parameters,'Enable','off');
%% PSF models
files=dir(fullfile(handles.smeagol_base_path,'components','PSF','*.m'));
models={'PSF model?'};
files([files.isdir]==0);
for k=1:length(files)
    [~,models{end+1}]=fileparts(files(k).name);
end    
set(handles.select_PSF_model,'String',models,'value',1);
set(handles.set_PSF_parameters,'Enable','off');
%% photophysics models
files=dir(fullfile(handles.smeagol_base_path,'components','photophys','*.m'));
models={'blink/bleach model?'};
files([files.isdir]==0);
for k=1:length(files)
    [~,models{end+1}]=fileparts(files(k).name);
end    
set(handles.select_photophysics_model,'String',models,'value',1);
set(handles.set_photophysics_parameters,'Enable','off');
%% background models
files=dir(fullfile(handles.smeagol_base_path,'components','background','*.m'));
models={'background?'};
files([files.isdir]==0);
for k=1:length(files)
    [~,models{end+1}]=fileparts(files(k).name);
end
set(handles.select_background_model,'String',models,'value',1);
set(handles.set_background_parameters,'Enable','off');
%% trj
% species information
set(handles.trj_numberOfSpecies,'Value',0);
set(handles.trj_numberOfSpecies,'String','0');
set(handles.table_species_diffusion,'Data',[]);
set(handles.trj_reactionFile,  'String','-');
set(handles.trj_trajectoryFile,'String','-');
%% output
set(handles.output_plotTrj,'Value',false);
set(handles.output_plotTrjZRange,'Enable','off');
set(handles.output_plotTrjZRange,'String','-1000 1000');
% reset ROI
set(handles.camera_pixLength_xRange_yRange,'Data',[0 0 0]');

% create options struct
%load('reset_opt','opt');
%handles.opt=opt;
handles.opt=struct('activation',struct,'baseIntensity',struct,'background',struct,'camera',struct,...
             'output',struct,'photophys',struct,'psf',struct,'sample',struct,...
             'trj',struct);
% save new GUI data
guidata(hObject, handles);
%% TRY TO ADD A BACKGROUND IMAGE
% This creates the 'background' axes
%ha = axes('units','normalized','position',[0 0 1 1]);
%mesord
% Move the background axes to the bottom
%uistack(ha,'bottom');
%
% Load in a background image and display it using the correct colors
% The image used below, is in the Image Processing Toolbox.  If you do not
% have %access to this toolbox, you can use another image file instead. 
%I=imread(fullfile(h.smeagol_base_path,'gui/Palantir_Stone_nocopyright.jpg'));
%imagesc(I);
%%colormap gray
%
% Turn the handlevisibility off so that we don't inadvertently plot into
% the axes again. Also, make the axes invisible.
% set(ha,'handlevisibility','off','visible','off')

% --- Outputs from this function are returned to the command line.
function varargout = SMeagol_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%uiwait(handles.figure1); % wait here until closing time

varargout{1} = handles.output;
%% internal functions
function opt=gui2options(hObject)
handles=guidata(hObject);
opt=struct;

% figure out the runinputfile path
[runinputroot,~,~]=fileparts(get(handles.str_runinputfile,'String'));

% trj
opt.trj.reactionFile    =   SM_relative_path_to_file(runinputroot,get(handles.trj_reactionFile,'String'));
opt.trj.trajectoryFile  =   SM_relative_path_to_file(runinputroot,get(handles.trj_trajectoryFile,'String'));
opt.trj.degradedName    =   get(handles.trj_degradedName,'String');
opt.trj.timeScale       =   str2num(get(handles.trj_timeScale,'String'));
opt.trj.voxelSize       =   str2num(get(handles.trj_voxelSize,'String'));

data=get(handles.table_species_diffusion,'Data');
N=size(data,1);
for k=1:N
    opt.trj.speciesNames{k}=data{k,1};
    opt.trj.D(k)=data{k,2};
end
set(handles.table_species_diffusion,'Data',data);

% activation
if(~isempty(handles.opt.activation) && ~isempty(handles.opt.activation.type) )
    opt.activation=handles.opt.activation;
else
    warndlg('No activation model chosen. Runinput file will be incomplete.')
end

% base intensity
if(~isempty(handles.opt.baseIntensity) && ~isempty(handles.opt.baseIntensity.type) )
    opt.baseIntensity=handles.opt.baseIntensity;
else
    warndlg('No base intensity model chosen. Runinput file will be incomplete.')
end

% PSF
if(~isempty(handles.opt.psf) && ~isempty(handles.opt.psf.type) )
    opt.psf=handles.opt.psf;
else
    warndlg('No psf model chosen. Runinput file will be incomplete.')
end

% photophysics
if(~isempty(handles.opt.photophys) && ~isempty(handles.opt.photophys.type) )
    opt.photophys=handles.opt.photophys;
else
    warndlg('No photophysics model chosen. Runinput file will be incomplete.')
end

% background
if(~isempty(handles.opt.background) && ~isempty(handles.opt.background.type) )
    opt.background=handles.opt.background;
else
    warndlg('No background model chosen. Runinput file will be incomplete.')
end

%  camera 
opt.camera.alpha  = str2num(get(handles.camera_alpha,'String'));
opt.camera.offset = str2num(get(handles.camera_offset,'String'));
opt.camera.sigmaReadout=str2num(get(handles.camera_sigmaReadout,'String'));

Data=get(handles.camera_pixLength_xRange_yRange,'Data');
opt.camera.pixLength=Data(1);
opt.camera.xrange_px=Data(2);
opt.camera.yrange_px=Data(3);

opt.camera.A=get(handles.camera_A,'Data');
opt.camera.b=get(handles.camera_b,'Data');

% output: [1x1 struct] : produce logical variables rather than doubles
opt.output.resultFile=SM_relative_path_to_file(runinputroot,get(handles.output_result_file,'String'));
opt.output.writeTifMovie    = get(handles.output_writeTifMovie,'Value')==1;
opt.output.plotTifMovie     = get(handles.output_plotTifMovie ,'Value')==1;
opt.output.showPhotons      = get(handles.output_showPhotons  ,'Value')==1;
opt.output.showEmitters     = get(handles.output_showEmitters ,'Value')==1;
opt.output.plotTrj          = get(handles.output_plotTrj      ,'Value')==1;

opt.output.movieLength      = round(str2double(get(handles.output_movieLength,'String')));
opt.output.maxFrames        = round(str2double(get(handles.output_maxFrames,'String')));

% defaults values waiting for implementation
opt.output.plotTrjZRange=str2num(get(handles.output_plotTrjZRange,'String'));

% true defaults:
opt.output.movieFormat='tiff';
opt.output.movieOptsImwrite={};

% sample: [1x1 struct]
opt.sample.dt=str2num(get(handles.sample_dt,'String'));
opt.sample.tE=str2num(get(handles.sample_tE,'String'));
function update_gui_settings(opt,handles)
% function to update all GUI parts according to a supplied options structure
% opt       - valid Palantir options structure
% handles   - list od GUI component handles

% trj
set(handles.trj_trajectoryFile,'String',fullfile(opt.runinputroot,opt.trj.trajectoryFile));
set(handles.trj_reactionFile,'String', fullfile(opt.runinputroot,opt.trj.reactionFile));
set(handles.trj_degradedName,'String', opt.trj.degradedName);
set(handles.trj_timeScale,'String',    num2str(opt.trj.timeScale));
set(handles.trj_voxelSize,'String',    num2str(opt.trj.voxelSize));

N=length(opt.trj.D);
set(handles.trj_numberOfSpecies,'String',int2str(N));
data=cell(N,2);
for k=1:N
    data{k,1}=opt.trj.speciesNames{k};
    data{k,2}=opt.trj.D(k);
end
set(handles.table_species_diffusion,'Data',data);

% output
set(handles.output_result_file,'String',fullfile(opt.runinputroot,opt.output.resultFile));
set(handles.output_writeTifMovie,'Value',opt.output.writeTifMovie);
set(handles.output_plotTifMovie, 'Value',opt.output.plotTifMovie);
set(handles.output_showPhotons,  'Value',opt.output.showPhotons);
set(handles.output_showEmitters, 'Value',opt.output.showEmitters);
set(handles.output_plotTrj,      'Value',opt.output.plotTrj);

set(handles.output_maxFrames,    'String',int2str(opt.output.maxFrames));
set(handles.output_movieLength,  'String',int2str(opt.output.movieLength));
set(handles.output_plotTrjZRange,'String',num2str(opt.output.plotTrjZRange));

% default values are used for these flags
%           movieFormat: 'tiff'
%      movieOptsImwrite: {}

% ROI/camera
set(handles.camera_alpha,'String',         num2str(opt.camera.alpha));
set(handles.camera_sigmaReadout,'String',  num2str(opt.camera.sigmaReadout));
set(handles.camera_offset,'String',        num2str(opt.camera.offset));

px_ROIx_ROIy=[opt.camera.pixLength;opt.camera.xrange_px;opt.camera.yrange_px];
set(handles.camera_pixLength_xRange_yRange,'Data',px_ROIx_ROIy);
set(handles.camera_A,'Data',opt.camera.A);
set(handles.camera_b,'Data',opt.camera.b);

% sampling
set(handles.sample_dt,'String',num2str(opt.sample.dt));
set(handles.sample_tE,'String',num2str(opt.sample.tE));

% copy activation settings if they are compatible with existing methods,
% complain otherwise.
actType=find(strcmpi(opt.activation.type,get(handles.select_activation_model,'String')));
if(isempty(actType))
    warndlg(['runinput file activity.type = ' opt.activation.type ' not recongized. Ignoring.'])
    set(handles.select_activation_model,'Value',1);
    set(handles.set_activation_parameters,'Enable','off');
    handles.opt.activation=struct('type','');
else
    set(handles.select_activation_model,'Value',actType);
    set(handles.set_activation_parameters,'Enable','on');
    handles.opt.activation=opt.activation;
end

% copy intensity settings if they are compatible with existing methods,
% complain otherwise. 
intType=find(strcmpi(opt.baseIntensity.type,get(handles.select_baseIntensity_model,'String')));
if(isempty(intType))
    warndlg(['runinput file baseIntensity.type = ' opt.baseIntensity.type ' not recongized.'])
    set(handles.select_baseIntensity_model,  'Value',1);
    set(handles.set_baseIntensity_parameters,'Enable','off');
    handles.opt.baseIntensity=struct('type','');
else
    set(handles.select_baseIntensity_model,  'Value',intType);
    set(handles.set_baseIntensity_parameters,'Enable','on');
    handles.opt.baseIntensity=opt.baseIntensity;
end

% copy PSF settings if they are compatible with existing methods, complain
% otherwise. 
psfType=find(strcmpi(opt.psf.type,get(handles.select_PSF_model,'String')));
if(isempty(psfType))
    warndlg(['runinput file psf.type = ' opt.psf.type ' not recongized.'])
    set(handles.select_PSF_model,'Value',1);
    set(handles.set_PSF_parameters,'Enable','off');
    handles.opt.psf=struct('type','');
else
    set(handles.select_PSF_model,'Value',psfType);
    set(handles.set_PSF_parameters,'Enable','on');
    handles.opt.psf=opt.psf;
end

% copy photophysics settings if they are compatible with existing methods, complain
% otherwise. 
photophysType=find(strcmpi(opt.photophys.type,get(handles.select_photophysics_model,'String')));
if(isempty(photophysType))
    warndlg(['runinput file photophys.type = ' opt.photophys.type ' not recongized.'])
    set(handles.select_photophysics_model,'Value',1);
    set(handles.set_photophysics_parameters,'Enable','off');
    handles.opt.photophys=struct('type','');
else
    set(handles.select_photophysics_model,'Value',photophysType);
    set(handles.set_photophysics_parameters,'Enable','on');
    handles.opt.photophys=opt.photophys;
end

% copy background settings if they are compatible with existing methods, complain
% otherwise. 
bgType=find(strcmpi(opt.background.type,get(handles.select_background_model,'String')));
if(isempty(bgType))
    warndlg(['runinput file psf.background = ' opt.background.type ' not recongized.'])
    set(handles.select_background_model,'Value',1);
    set(handles.set_background_parameters,'Enable','off');
    handles.opt.background=struct('type','');
else
    set(handles.select_background_model,'Value',bgType);
    set(handles.set_background_parameters,'Enable','on');
    handles.opt.background=opt.background;
end

guidata(handles.figure1,handles) % store updated values
%% callback functions for flourescence parameters
function select_activation_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_activation_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns select_activation_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_activation_model

% selected activation model
contents = cellstr(get(hObject,'String'));
selectedActivationModel=contents{get(hObject,'Value')};
% activate parameters if file exists, deactivate otherwise
if(exist(selectedActivationModel,'file')~=2)
    set(handles.set_activation_parameters,'Enable','off')    
    newActivationOpt=struct('type','');
else    
    set(handles.set_activation_parameters,'Enable','on')
    newActivationOpt=eval(selectedActivationModel); % get default parameters for new model
    newActivationOpt.type=selectedActivationModel;  % set name of activation model
end
% if the selection is new, discard old parameters
oldActivationOpt=handles.opt.activation;
if(isfield(oldActivationOpt,'type') && strcmp(oldActivationOpt.type,newActivationOpt.type) )
    % then selection was not changed, and options should not be updated
else
    handles.opt.activation=newActivationOpt;
end
disp('current activation options:')
handles.opt.activation
guidata(hObject,handles); % store changes in handles 
function select_baseIntensity_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_baseIntensity_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns select_baseIntensity_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_baseIntensity_model

% selected activation model
contents = cellstr(get(hObject,'String'));
selectedBaseIntensityModel=contents{get(hObject,'Value')};
% activate parameters if file exists, deactivate otherwise
if(exist(selectedBaseIntensityModel,'file')~=2)
    set(handles.set_baseIntensity_parameters,'Enable','off')    
    newBaseIntensityOpt=struct('type','');
else    
    set(handles.set_baseIntensity_parameters,'Enable','on')
    newBaseIntensityOpt=eval(selectedBaseIntensityModel); % get default parameters for new model
    newBaseIntensityOpt.type=selectedBaseIntensityModel;  % set name of activation model
end
% if the selection is new, discard old parameters
oldBaseIntensityOpt=handles.opt.baseIntensity;
if(isfield(oldBaseIntensityOpt,'type') && strcmp(oldBaseIntensityOpt.type,newBaseIntensityOpt.type) )
    % then selection was not changed, and options should not be updated
else
    handles.opt.baseIntensity=newBaseIntensityOpt;
end
disp('current base intensity options:')
handles.opt.baseIntensity
guidata(hObject,handles); % store changes in handles 
function select_PSF_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_PSF_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns select_PSF_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_PSF_model

% selected PSF model
contents = cellstr(get(hObject,'String'));
selectedPSFmodel=contents{get(hObject,'Value')};
% activate parameters if file exists, deactivate otherwise
if(exist(selectedPSFmodel,'file')~=2)
    set(handles.set_PSF_parameters,'Enable','off')    
    newPSFopt=struct('type','');
else    
    set(handles.set_PSF_parameters,'Enable','on')
    newPSFopt=eval(selectedPSFmodel);
    newPSFopt.type=selectedPSFmodel;
end
% if the selection is new, discard old parameters
oldPSFopt=handles.opt.psf;
if(isfield(oldPSFopt,'type') && strcmp(oldPSFopt.type,newPSFopt.type) )
    % then selection was not changed, and options should not be updated
else
    handles.opt.psf=newPSFopt;
end
disp('current PSF options:')
handles.opt.psf
guidata(hObject,handles); % store changes in handles 
function select_photophysics_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_photophysics_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns select_photophysics_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_photophysics_model

% selected background (_photophys_) model
contents = cellstr(get(hObject,'String'));
selected_photophys_model=contents{get(hObject,'Value')};
% activate parameters if file exists, deactivate otherwise
if(exist(selected_photophys_model,'file')~=2)
    set(handles.set_photophysics_parameters,'Enable','off')    
    new_photophys_opt=struct('type','');
else    
    set(handles.set_photophysics_parameters,'Enable','on')
    new_photophys_opt=eval(selected_photophys_model);
    new_photophys_opt.type=selected_photophys_model;
end
% if the selection is new, discard old parameters
old_photophys_opt=handles.opt.photophys;
if(isfield(old_photophys_opt,'type') && strcmp(old_photophys_opt.type,new_photophys_opt.type) )
    % then selection was not changed, and options should not be updated
else
    handles.opt.photophys=new_photophys_opt;
end
disp('current photophysics options:')
handles.opt.photophys
guidata(hObject,handles); % store changes in handles 
function select_background_model_Callback(hObject, eventdata, handles)
% hObject    handle to select_background_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns select_background_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from select_background_model

% selected background (_photophys_) model
contents = cellstr(get(hObject,'String'));
selected_photophys_model=contents{get(hObject,'Value')};
% activate parameters if file exists, deactivate otherwise
if(exist(selected_photophys_model,'file')~=2)
    set(handles.set_background_parameters,'Enable','off')    
    new_photophys_opt=struct('type','');
else    
    set(handles.set_background_parameters,'Enable','on')
    new_photophys_opt=eval(selected_photophys_model);
    new_photophys_opt.type=selected_photophys_model;
end
% if the selection is new, discard old parameters
old_photophys_opt=handles.opt.background;
if(isfield(old_photophys_opt,'type') && strcmp(old_photophys_opt.type,new_photophys_opt.type) )
    % then selection was not changed, and options should not be updated
else
    handles.opt.background=new_photophys_opt;
end
disp('current background options:')
handles.opt.background
guidata(hObject,handles); % store changes in handles 
function set_activation_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to set_activation_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cellString=get(handles.select_activation_model,'String');
actType=cellString{get(handles.select_activation_model,'Value')};
try % we do not want the user to see or change the 'type' field
    handles.opt.activation=rmfield(handles.opt.activation,'type');
catch
end
newACT=component_gui3(actType,handles.opt.activation);
newACT.type=actType;        % always save what type of PSF model is to be used
handles.opt.activation=newACT;     % save options in GUI data structure
disp('new activation parameters')
disp(handles.opt.activation)
guidata(hObject,handles);   % store changes in handles 
function set_baseIntensity_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to set_baseIntensity_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellString=get(handles.select_baseIntensity_model,'String');
bimType=cellString{get(handles.select_baseIntensity_model,'Value')};
try % we do not want the user to see or change the 'type' field
    handles.opt.baseIntensity=rmfield(handles.opt.baseIntensity,'type');
catch
end
newBIM=component_gui3(bimType,handles.opt.baseIntensity);
newBIM.type=bimType;        % always save what type of PSF model is to be used
handles.opt.baseIntensity=newBIM;     % save options in GUI data structure
disp('new base intensity parameters')
disp(handles.opt.baseIntensity)
guidata(hObject,handles);   % store changes in handles 
function set_PSF_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to set_PSF_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellString=get(handles.select_PSF_model,'String');
psfType=cellString{get(handles.select_PSF_model,'Value')};
try % we do not want the user to see or change the 'type' field
    handles.opt.psf=rmfield(handles.opt.psf,'type');
catch
end
newPSF=component_gui3(psfType,handles.opt.psf);
newPSF.type=psfType;        % always save what type of PSF model is to be used
handles.opt.psf=newPSF;     % save options in GUI data structure
disp('new PSF parameters')
disp(handles.opt.psf)
guidata(hObject,handles);   % store changes in handles 
function set_background_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to set_background_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellString=get(handles.select_background_model,'String');
backgroundType=cellString{get(handles.select_background_model,'Value')};
try % we do not want the user to see or change the 'type' field
    handles.opt.background=rmfield(handles.opt.background,'type');
catch
end
newBackground=component_gui3(backgroundType,handles.opt.background);
newBackground.type=backgroundType;        % always save what type of PSF model is to be used
handles.opt.background=newBackground;     % save options in GUI data structure
disp('new background parameters')
disp(handles.opt.background)
guidata(hObject,handles);   % store changes in handles 
function set_photophysics_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to set_photophysics_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cellString=get(handles.select_photophysics_model,'String');
photophysType=cellString{get(handles.select_photophysics_model,'Value')};
try % we do not want the user to see or change the 'type' field
    handles.opt.photophys=rmfield(handles.opt.photophys,'type');
catch
end
newPhotophys=component_gui3(photophysType,handles.opt.photophys);
newPhotophys.type=photophysType;        % always save what type of PSF model is to be used
handles.opt.photophys=newPhotophys;     % save options in GUI data structure
disp('new photophysics parameters')
disp(handles.opt.photophys)
guidata(hObject,handles);   % store changes in handles 
% --- Executes during object creation, after setting all properties.
function select_activation_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_activation_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function select_baseIntensity_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_baseIntensity_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function select_PSF_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_PSF_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function select_photophysics_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_photophysics_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function select_background_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to select_background_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% callback functions for camera ROI parameters
function camera_set_A_b_gui_Callback(hObject, eventdata, handles)
% hObject    handle to camera_set_A_b_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% check for trj file name
trjFileName=get(handles.trj_trajectoryFile,'String');
if( ~exist(trjFileName,'file'))
    errordlg('Trajectory file needed to start the A,b GUI. Cannot find it.')
    return;
end

% check if a GUI already exists, and ask to delete it
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
    deleteOldABgui=questdlg('An A,b GUI seems to be already in existence. Close it?','','No');
    if(strcmp(deleteOldABgui,'Yes'))
        % try to delete the old GUI and its figure window
        try
            ABgui.closeThisGUI();
            clear ABgui
        catch % if something went wrong , then perhaps it was no functional GUI?            
        end
        set(handles.camera_set_A_b_gui,'UserData',[]); % delete the handle
    else
        return;
    end
end

% call ABgui and save a handle to it
ABgui=palantirAB_gui(handles.figure1);
set(handles.camera_set_A_b_gui,'UserData',ABgui);
guidata(hObject,handles);   % store changes in handles 
function camera_pixLength_xRange_yRange_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to camera_pixLength_xRange_yRange (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%%% should not need to get new handles?
%handles=guidata(hObject);
% if an AB_gui exists, force an update of it
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
   ABgui.updateROI();
end 
%%%guidata(hObject,handles);   % store changes in handles 
function camera_A_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to camera_A (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% if an AB_gui exists, force an update of it
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
   ABgui.updateROI();
end    
%%%guidata(hObject,handles);   % store changes in handles 
function camera_b_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to camera_b (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% if an AB_gui exists, force an update of it
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
   ABgui.updateROI();
end
%%%guidata(hObject,handles);   % store changes in handles 
%% callback functions for camera noise parameters
function camera_sigmaReadout_Callback(hObject, eventdata, handles)
% hObject    handle to camera_sigmaReadout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_sigmaReadout,OK]=str2num(get(handles.camera_sigmaReadout,'String'));
if(OK && length(new_sigmaReadout)==1 && new_sigmaReadout>=0)    
    set(handles.camera_sigmaReadout,'String',new_sigmaReadout);
else
    warndlg('The EMCCD readout noise std should be a non-negative number.')
    set(handles.camera_sigmaReadout,'String',num2str(handles.opt.camera.sigmaReadout));
end
guidata(hObject,handles); % save new values
function camera_offset_Callback(hObject, eventdata, handles)
% hObject    handle to text_camera_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_offset,OK]=str2num(get(handles.camera_offset,'String'));
new_offset=round(new_offset);
if(OK && length(new_offset)==1 && new_offset>=0)
    set(handles.camera_offset,'String',int2str(new_offset));
else
    warndlg('The EMCCD offset should be a non-negative integer.')
    set(handles.camera_offset,'String','-');
end
guidata(hObject,handles); % save new values
function camera_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to camera_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_alpha,OK]=str2num(get(handles.camera_alpha,'String'));
if(OK && length(new_alpha)==1 && new_alpha>0)
    set(handles.camera_alpha,'String',num2str(new_alpha));
else
    set(handles.camera_alpha,'String','-');
    warndlg('The EMCCD inverse gain must be a positive number.')
end
guidata(hObject,handles); % save new value
% --- Executes during object creation, after setting all properties.
function camera_sigmaReadout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_sigmaReadout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function camera_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function camera_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to camera_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
%
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function text_camera_offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_camera_offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% callbacks for illumination/exposure
function sample_dt_Callback(hObject, eventdata, handles)
% hObject    handle to sample_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_dt as text
%        str2double(get(hObject,'String')) returns contents of sample_dt as a double
[new_sample_dt,OK]=str2num(get(handles.sample_dt,'String'));
if(OK && length(new_sample_dt)==1 && new_sample_dt>0)
    set(handles.sample_dt,'String',num2str(new_sample_dt));
else
    warndlg('The sampling time should be a positive number.')
    set(handles.sample_dt,'String','-');
    return;
end
[old_sample_tE,OK]=str2num(get(handles.sample_tE,'String'));
if(OK && new_sample_dt < old_sample_tE)
    warndlg('Setting aquisition time = sampling time for consistency.')
    set(handles.sample_tE,'String',get(handles.sample_dt,'String'));
end
guidata(hObject,handles); % save new values
function sample_tE_Callback(hObject, eventdata, handles)
% hObject    handle to sample_tE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_sample_tE,OK]=str2num(get(handles.sample_tE,'String'));
if(OK && length(new_sample_tE)==1 && new_sample_tE>0)
    set(handles.sample_tE,'String',num2str(new_sample_tE));
else
    warndlg('The aquisition time should be a positive number.')
    set(handles.sample_tE,'String','-');
    return
end
[old_sample_dt,OK]=str2num(get(handles.sample_dt,'String'));
if(OK && old_sample_dt < new_sample_tE)
    warndlg('Setting sampling time = aquisition time for consistency.')
    set(handles.sample_dt,'String',get(handles.sample_tE,'String'));
end
guidata(hObject,handles); % save new values
% --- Executes during object creation, after setting all properties.
function sample_dt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_dt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function sample_tE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_tE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% callback functions for input data
function select_trj_file_Callback(hObject, eventdata, handles)
% hObject    handle to select_trj_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[trjpath,~,~]=fileparts(get(handles.trj_trajectoryFile,'String'));
[rctpath,~,~]=fileparts(get(  handles.trj_reactionFile,'String'));
if(    ~isempty(trjpath) && exist(trjpath,'dir'))
    cd(trjpath)
elseif(~isempty(rctpath) && exist(rctpath,'dir'))
    cd(rctpath)
end
[trjfile,trjpath]=uigetfile({'*'},'select input trajectory file');
cd(handles.gui_start_path)
set(handles.trj_trajectoryFile,'String',fullfile(trjpath,trjfile));

% if an AB_gui exists, force an update of data
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))   
   %ABgui.updateROI();
   %%% needs testing!
   ABgui.updateTrj();   
end
% keep GUIdata handles up to date
guidata(hObject,handles)
function select_rct_file_Callback(hObject, eventdata, handles)
% hObject    handle to select_rct_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[trjpath,~,~]=fileparts(get(handles.trj_trajectoryFile,'String'));
[rctpath,~,~]=fileparts(get(handles.trj_reactionFile,  'String'));
if(    ~isempty(trjpath) && exist(trjpath,'dir'))
    cd(trjpath)
elseif(~isempty(rctpath) && exist(rctpath,'dir'))
    cd(rctpath)
end
[rctfile,rctpath]=uigetfile({'*'},'select input reaction file');
cd(handles.gui_start_path)
set(handles.trj_reactionFile,'String',fullfile(rctpath,rctfile));
guidata(hObject,handles); % save new values
function read_mesoRDlog_Callback(hObject, eventdata, handles)
% hObject    handle to read_mesoRDlog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% start in directory of trajectory or reaction files
[trjpath,~,~]=fileparts(get(handles.trj_trajectoryFile,'String'));
[rctpath,~,~]=fileparts(get(  handles.trj_reactionFile,'String'));
if(    ~isempty(trjpath) && exist(trjpath,'dir'))
    cd(trjpath)
elseif(~isempty(rctpath) && exist(rctpath,'dir'))
    cd(rctpath)
end
[logfile,logpath]=uigetfile(...
    {'molecule_tracking_log.txt','mesoRD log (molecule_tracking_log.txt)';
    '*.txt','text file (.txt)';'*','All files'},'select mesoRD log file');
cd(handles.gui_start_path);
if(~ischar(logfile)) % failed or cancelled logfile selection
    return
end
guidata(hObject,handles); % save new values

% parse the log file
readFromMesoRDlogFile(hObject,fullfile(logpath,logfile));
function readFromMesoRDlogFile(hObject,logPathFile)
[trj,exitflag]=SM_parse_mesoRD_logfile(logPathFile);
handles=guidata(hObject);
if(exitflag==0)
   %set(handles.trj_timestep,'String',num2str(trj.timestep));
   set(handles.trj_timeScale,'String',num2str(trj.timeScale));
   set(handles.trj_voxelSize,'String',num2str(trj.voxelSize));
   
   N=length(trj.D);
   set(handles.trj_numberOfSpecies,'String',int2str(N));
   
   data=cell(N,2);
   for k=1:N
      data{k,1}=trj.speciesNames{k};
      data{k,2}=trj.D(k);
   end
   set(handles.table_species_diffusion,'Data',data);
end
% possibly update reaction and trajectory files as well?
b=questdlg('Look for matching trajectory and reaction files?','Yes','No');
[logpath,logfile,~]=fileparts(logPathFile);
if(strcmp(b,'Yes'))
   new_trjFile=fullfile(logpath,'trajectories.txt');
   new_rctFile=fullfile(logpath,'reactions.txt');
   if( exist(new_trjFile,'file') && exist(new_rctFile,'file') )
       set(handles.trj_reactionFile  ,'String',new_rctFile);
       set(handles.trj_trajectoryFile,'String',new_trjFile);
       % if an AB_gui exists, force an update of data
       ABgui=get(handles.camera_set_A_b_gui,'UserData');
       if(~isempty(ABgui))
           ABgui.updateTrj();
           %ABgui.updateROI();
       end       
   else
       warndlg('Could not find trajectories.txt and reactions.txt in the lof file folder.')
   end
       
end 
guidata(hObject,handles);
function trj_timeScale_Callback(hObject, eventdata, handles)
% hObject    handle to trj_timeScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trj_timeScale as text
%        str2double(get(hObject,'String')) returns contents of trj_timeScale as a double
newTimeScale=str2num(get(hObject,'String'));
if( isnan(newTimeScale) || newTimeScale <=0 || ~isreal(newTimeScale))
    errordlg('The time scale factor must be a positive real number.')
    set(hObject,'String','-')
else
    set(hObject,'String',num2str(newTimeScale));
end
%guidata(hObject,handles);   % store changes in handles 
function trj_voxelSize_Callback(hObject, eventdata, handles)
% hObject    handle to trj_voxelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trj_voxelSize as text
%        str2double(get(hObject,'String')) returns contents of trj_voxelSize as a double
newVoxelSize=str2num(get(hObject,'String'));
if( isnan(newVoxelSize) || newVoxelSize <=0 || ~isreal(newVoxelSize))
    errordlg('Simulation voxel size must be a positive real number')
    set(hObject,'String','-')
else
    set(hObject,'String',num2str(newVoxelSize));
end
% if an AB_gui exists, force an update of it
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
   ABgui.updateROI();
end
function edit_species_list_Callback(hObject, eventdata, handles)
% hObject    handle to edit_species_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function trj_degradedName_Callback(hObject, eventdata, handles)
% hObject    handle to trj_degradedName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
new_trj_degradedName=get(hObject,'String');
% remove spaces 
ind=find(new_trj_degradedName~=' ');
if(length(ind)~=length(new_trj_degradedName) || isempty(ind))
    warndlg('Degraded species name must be a nonempty string without spaces.')
    new_trj_degradedName=new_trj_degradedName(new_trj_degradedName~=' ');
    if(isempty(new_trj_degradedName))
        new_trj_degradedName='-1';
    end
    set(hObject,'String',new_trj_degradedName);    
end
function trj_numberOfSpecies_Callback(hObject, eventdata, handles)
% hObject    handle to trj_numberOfSpecies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trj_numberOfSpecies as text
%        str2double(get(hObject,'String')) returns contents of trj_numberOfSpecies as a double
%[oldN,~]=size(handles
[oldN,~]=size(get(handles.table_species_diffusion,'Data'));
newN=round(str2num(get(hObject,'String')));
if( ~isempty(newN) && isfinite(newN) && isreal(newN) && newN>=0 )
    set(hObject,'String',int2str(newN));
    data=get(handles.table_species_diffusion,'Data');    
    if(newN<oldN)
        data=data(1:newN,:);        
    elseif(newN>oldN)        
        for k=oldN+1:newN
            data{k,1}='-';
            data{k,2}=0;
        end
    end
    set(handles.table_species_diffusion,'Data',data);
else
    errordlg('Need a non-negative integer number of species')
end
guidata(hObject,handles);   % store changes in handles 
% --- Executes during object creation, after setting all properties.
function trj_degradedName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trj_degradedName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function trj_timeScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trj_timeScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function trj_voxelSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trj_voxelSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function trj_numberOfSpecies_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trj_numberOfSpecies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% callback functions for runinput file, GUI, simulation
function set_runinput_file_Callback(hObject, eventdata, handles)
% hObject    handle to set_runinput_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% suggest that runinput file is colocalized with output
%[respath,~,~]=fileparts(  handles.output_result_file.String);

% better: start in current matlab folder
[respath,~,~]=fileparts(fullfile(pwd,'.'));

if(    ~isempty(respath) && exist(respath,'dir'))
    cd(respath)
end
[filename, pathname] = uiputfile('*.m','select runinput file to save to');
cd(handles.gui_start_path)

if(isnumeric(filename)) % then file selection was aborted/broken/or cancelled
    return;
end
% if not, then we have a runinput file candidate
set(handles.str_runinputfile,'String',fullfile(pathname,filename));
set(handles.save_runinput_file,'Enable','on');
set(handles.save_and_run,'Enable','on');
guidata(hObject,handles);   % store changes in handles 
function save_runinput_file_Callback(hObject, eventdata, handles)
% hObject    handle to save_runinput_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% extract path to runinput file
runinputfile=get(handles.str_runinputfile,'String');

% extract options
opt=gui2options(hObject);

% change explicit paths to relative paths ?
% ML: not sure how to do that.

flag=SM_writeRuninputFile(opt,runinputfile);
if(flag>=0)
   disp(['Save runinput file in ' runinputfile ])
else
    warndlg('Something went wrong during runiput file creation.')
    %%% TBA: save the state of the GUI for debugging
end
function save_and_run_Callback(hObject, eventdata, handles)
% hObject    handle to save_runinput_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
save_runinput_file_Callback(hObject, eventdata, handles);

% run simulation, using the current options% extract path to runinput file
runinputfile=get(handles.str_runinputfile,'String');

% extract options

% run simulation, using the current options: somehow, Matlab has
% difficulties reading the latest updated runinputfile, and rehash does not
% help
opt=gui2options(hObject);
opt2=SM_getOptions(runinputfile);
opt.runinputfile=opt2.runinputfile;
opt.runinputroot=opt2.runinputroot;
SM_runsimulation(opt);

function load_runinput_parameters_Callback(hObject, eventdata, handles)
% hObject    handle to load_runinput_parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('load new runinput parameters from file')
[optfile,optpath]=uigetfile('*.m','Select a Palantir runinput file.');
optfile=fullfile(optpath,optfile);
if(exist(optfile,'file'))    
    % read options
    try
        opt=SM_getOptions(optfile);
    catch
        errordlg(['Could not extract valid Palantir options from ' ...
            optfile ]);
        return
    end
    % some validation?                
    
    % if an AB_gui exists, force an update of data
    ABgui=get(handles.camera_set_A_b_gui,'UserData');
    if(~isempty(ABgui))
        ABgui.updateTrj();
        %ABgui.updateROI();
    end
end
% update GUI state with new options
update_gui_settings(opt,handles);
% no further guidata udate needed, that would overwrite what
% update_gui_setting just did
function close_gui_Callback(hObject, eventdata, handles)
% hObject    handle to close_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ABgui=get(handles.camera_set_A_b_gui,'UserData');
if(~isempty(ABgui))
    % close external AB-gui if it exists
    ABgui.closeThisGUI();
end
delete(handles.figure1); % continue execute the outputFcn
% --- Executes during object creation, after setting all properties.
%% callback for output fields
function output_movieLength_Callback(hObject, eventdata, handles)
% hObject    handle to output_movieLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_output_movieLength,OK]=str2num(get(handles.output_movieLength,'String'));
if(OK && length(new_output_movieLength)==1 && new_output_movieLength>0)
    set(handles.output_movieLength,'String',int2str(round(new_output_movieLength)));
else
    warndlg('The movieLength should be a positive integer.')
    set(handles.output_movieLength,'String','-');
end
guidata(hObject,handles); % save new values
function output_maxFrames_Callback(hObject, eventdata, handles)
% hObject    handle to output_maxFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_output_maxFrames,OK]=str2num(get(handles.output_maxFrames,'String'));
if(OK && length(new_output_maxFrames)==1 && new_output_maxFrames>0)
    set(handles.output_maxFrames,'String',int2str(round(new_output_maxFrames)));
else
    warndlg('The maxFrames should be a positive integer.')
    set(handles.output_maxFrames','String','-');
end
guidata(hObject,handles); % save new values
function set_result_file_Callback(hObject, eventdata, handles)
% hObject    handle to set_result_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% suggest that runinput file is colocalized with output
[RIpath,~,~]=fileparts(get(handles.str_runinputfile,'String'));
if(    ~isempty(RIpath) && exist(RIpath,'dir'))
    cd(RIpath)
end
[filename, pathname] = uiputfile('*.mat','select runinput file to save to');
cd(handles.gui_start_path)

if(isnumeric(filename)) % then file selection was aborted/broken/or cancelled
    return;
end
set(handles.output_result_file,'String',fullfile(pathname,filename));
function output_plotTrjZRange_Callback(hObject, eventdata, handles)
% hObject    handle to output_plotTrjZRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[new_output_plotTrjZRange,OK]=str2num(get(handles.output_plotTrjZRange,'String')); %#ok<ST2NM>
if(OK && length(new_output_plotTrjZRange)==2 && diff(new_output_plotTrjZRange)>0)    
        set(handles.output_plotTrjZRange,'String',num2str(new_output_plotTrjZRange));
    else
    warndlg('The plotTrjZRange should be a 2-vector [zmin zmax], with zmin < zmax.')
    set(handles.output_plotTrjZRange,'String','-1 1');
end
guidata(hObject,handles); % save new values
function output_plotTrj_Callback(hObject, eventdata, handles)
h=guidata(hObject);
if(get(hObject,'Value'))
   set(h.output_plotTrjZRange,'Enable','on');
else
   set(h.output_plotTrjZRange,'Enable','off');
end
guidata(hObject,h);   % store changes in handles 
function output_writeTifMovie_Callback(hObject, eventdata, handles)
function output_plotTifMovie_Callback(hObject, eventdata, handles)
function output_showPhotons_Callback(hObject, eventdata, handles)
function output_showEmitters_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function output_plotTrjZRange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_plotTrjZRange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_movieLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_movieLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function output_maxFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output_movieLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run_mesoRD_button.
function run_mesoRD_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_mesoRD_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%% need to make things happen here!!!
h=guidata(hObject);
logfile=SM_mesoRD_gui;

if(~exist(logfile,'file'))
    return
end
b=questdlg('Read trj info from target logfile?','Yes','No');
if(strcmp(b,'Yes'))
    readFromMesoRDlogFile(hObject,logfile);
end

% --- Executes during object creation, after setting all properties.
function run_mesoRD_button_CreateFcn(hObject, eventdata, handles)
% hObject    handle to run_mesoRD_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
