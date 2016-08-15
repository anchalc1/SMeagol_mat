function varargout = palantirAB_gui(varargin)
%PALANTIRAB_GUI M-file for palantirAB_gui.fig
%      PALANTIRAB_GUI, by itself, creates a new PALANTIRAB_GUI or raises the existing
%      singleton*.
%
%      H = PALANTIRAB_GUI returns the handle to a new PALANTIRAB_GUI or the handle to
%      the existing singleton*.
%
%      PALANTIRAB_GUI('Property','Value',...) creates a new PALANTIRAB_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to palantirAB_gui_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PALANTIRAB_GUI('CALLBACK') and PALANTIRAB_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PALANTIRAB_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help palantirAB_gui

% Last Modified by GUIDE v2.5 13-Jan-2015 15:38:52

%% copyright notice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% palantirAB_gui.m, sub-GUI in the SMeagol package
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
                   'gui_OpeningFcn', @palantirAB_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @palantirAB_gui_OutputFcn, ...
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

% --- Executes just before palantirAB_gui is made visible.
function palantirAB_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% UIWAIT makes palantirAB_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%% set some default values
set(handles.str_angle_increment,'String','30')
set(handles.str_angle_increment,'Value',30)
set(handles.data_trjpts_stride,'Data',[500 10]); % 500 points, every 10 pts

% save handle to the parent GUI
handles.parentGUI_hObject=varargin{1};
%% set non-interruptible to decrease mess...
set(handles.figure1,'Interruptible','off');
%% open and initialize plots
x=[0 0];y=[0 0];z=[0 0];xL=1;yL=1; % initial values, to be updated later

handles.plotWindow=struct;
handles.plotWindow.figure=figure;
handles.plotWindow.xTrj=[x;y;z];
set(handles.plotWindow.figure,'CloseRequestFcn',@(s,c)(palantirAB_gui_CloseFcn(s,c,hObject)));
set(handles.figure1,'CloseRequestFcn',@(s,c)(palantirAB_gui_CloseFcn(s,c,hObject)));

clf
% xy-plane with red ROI
handles.plotWindow.hxy=subplot(2,2,1);
handles.plotWindow.dat=plot(x,y,'k.');
hold on
handles.plotWindow.ROIpatch=patch(xL*[0 1 1 0 0],yL*[0 0 1 1 0],'r','edgecol','none','facealpha',0.5);
xlabel('x'),ylabel('y')
axis equal
grid on

% xz-plane
handles.plotWindow.hxz=subplot(2,2,3);
handles.plotWindow.dat(2)=plot(x,z,'k.');
hold on
xlabel('x'),ylabel('z')
axis equal
grid on
% xz-plane
handles.plotWindow.hyz=subplot(2,2,4);
handles.plotWindow.dat(3)=plot(y,z,'k.');
hold on
xlabel('y'),ylabel('z')
axis equal
grid on

% add some intstructions in the plot window
subplot(2,2,2)
axis([0 1 0 1])
axis off

msg={'The graphs show the present',...
     'ROI in three views. The red',...
     'rectangle is the region to',...
     'be imaged (ROI), z=0 is the',...
     'simulated focal plane, black',...
     'dots show some positions',...
     'from the underlying trajectory,',...
     'indicating its geometry.'};

text(0,1,msg,'VerticalAlignment','top');
%% add moveable rectangle in all three planes
ROIrect(1)=imrect(handles.plotWindow.hxy,[0 0 xL yL]);
ROIrect(2)=imrect(handles.plotWindow.hxz,[0 0 xL 0]);
ROIrect(3)=imrect(handles.plotWindow.hyz,[0 0 yL 0]);
% fixed size
setResizable(ROIrect(1),false)
setResizable(ROIrect(2),false)
setResizable(ROIrect(3),false)
handles.plotWindow.ROIrect=ROIrect;

% callback functions for moving them
%addNewPositionCallback(ROIrect(1),@(pos)(updateROIxy(pos,hObject)));
%addNewPositionCallback(ROIrect(2),@(pos)(updateROIxz(pos,guidata(hObject))));
%addNewPositionCallback(ROIrect(3),@(pos)(updateROIyz(pos,guidata(hObject))));
%% Update handles structure, ROI, and trajectory data
guidata(hObject, handles);
% read data 
updateTrjFileAndReplot(hObject);
% --- Outputs from this function are returned to the command line.
function varargout = palantirAB_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
outStruct=struct;
outStruct.closeThisGUI=@()(palantirAB_gui_CloseFcn([],[],hObject));
outStruct.updateROI=@()(updateROI_fromMainGUI(hObject));
outStruct.updateTrj=@()(updateTrjFileAndReplot(hObject));

varargout{1} = outStruct;
%% internal functions
function palantirAB_gui_CloseFcn(src,callbackdata,hObject)
% Close request function
% to display a question dialog box
try % try to close both the GUI and the ROI plot window
    handles=guidata(hObject);
    delete(handles.plotWindow.figure);
    delete(handles.figure1)
    fclose(handles.plotWindow.trjFid);
catch % something is wrong (perhaps some figure is already deleted or not opened) 
      % simply close the current window
    delete(gcf)
    return
end


% and also delete the parent GUI handle to the A,B gui (= this GUI).
try
    ph=guidata(handles.parentGUI_hObject);
    set(ph.camera_set_A_b_gui,'UserData',[]);
catch
    % if this does not work, just ignore it (the main GUI might be deleted
    % already).
end
function updateTrjFileAndReplot(hObject)
handles=guidata(hObject);
% close existing trajectory file if not done yet
if(isfield(handles.plotWindow,'trjFid'))
    try
        fclose(handles.plotWindow.trjFid)
    catch
    end
end
% reset old information before trying to add new
handles.plotWindow.trjFid=[];
handles.plotWindow.trjSpcStr=[];

% now open the current file
ph=guidata(handles.parentGUI_hObject);
fid=fopen(get(ph.trj_trajectoryFile,'String'),'r');
Data=get(ph.table_species_diffusion,'Data'); % names of the species
spcStr=Data(:,1); % names of the species
if(fid==-1)
    errordlg(['ABgui could not open the trjectory file: ' ph.trj_trajectoryFile])
    return
end
[~,~,~,x]=SM_readLine(fid,eye(3),zeros(3,1),1,spcStr);
if(isempty(x)) % check for readable positions, then close and reopen
    errordlg(['ABgui found empty trj file: ' ph.trj_trajectoryFile])
    return
end
fclose(fid); % close and reopen 
fid=fopen(get(ph.trj_trajectoryFile,'String'),'r');

% save file id and species str.
handles.plotWindow.trjFid=fid;
handles.plotWindow.spcStr=spcStr;
guidata(handles.figure1,handles);
% update the displayed data
updateTrjData(hObject)
function updateTrjData(hObject)
% replace the displayed trajectory data with new positions from the old
% trajectory file, starting from the beginning if necessary.

hh=msgbox('Reading trajectory data (this message will self-delete when done).');
handles=guidata(hObject);
fid=handles.plotWindow.trjFid;
spcStr=handles.plotWindow.spcStr;
Data=get(handles.data_trjpts_stride,'Data');
numPtsMax=Data(1);
Pacc=1/Data(2); % position retention probability

X=[];
numPts=0;
fileLoops=0;
ph=guidata(handles.parentGUI_hObject);
while(numPts <= numPtsMax)
       [~,~,~,x]=SM_readLine(fid,eye(3),zeros(3,1),1,spcStr);
       if(isempty(x)) % then close and reread from the beginning
           fclose(fid);
           fid=fopen(get(ph.trj_trajectoryFile,'String'),'r');
           [~,~,~,x]=SM_readLine(fid,eye(3),zeros(3,1),1,spcStr);
           fileLoops=fileLoops+1;
       end
       if(isempty(X) && fileLoops > 1) % then close and reread from the beginning
           errordlg('trj file appears to be empty.')
           X=zeros(3,2);
           break
       end
       ind=find(rand(size(x,1),1)<Pacc);
       X=[X x(ind,:)'];
       numPts=numPts+length(ind);
end
handles.plotWindow.xTrj=X(:,1:numPtsMax);
delete(hh);
% updated handles fields
guidata(handles.figure1,handles);
updateROI_fromMainGUI(handles.figure1);
function updateROI_fromMainGUI(hObject)
% read ROI parameters from the main GUI and update the various plots.
% This function is thought to be called when ROI-related properties are
% changed in the main GUI. 

% update this GUI data (in case this function is called with outdated
% handles, for example by the inline construction sent to the main GUI
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

% ROI parameters from the parent GUI
vxl=str2double(get(ph.trj_voxelSize,'String'));
camA=get(ph.camera_A,'Data');
camB=get(ph.camera_b,'Data');
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
pixLength= Data(1);
xLpx     = Data(2);
yLpx     = Data(3);

% ROI: Images are defined with positions in the middle of pixels. So,
% x=1,y=1 should be in the middle of the (1,1) pixel, hence the +0.5
xPatch=pixLength*(+0.5+[0 xLpx xLpx 0 0]);
yPatch=pixLength*(+0.5+[0 0 yLpx yLpx 0]);
zPatch=[0 0 0 0 0];
set(handles.plotWindow.ROIpatch,'xdata',xPatch,'ydata',yPatch);
% update rectangle properties (xy,xz,yz)
handles.plotWindow.ROIrect(1).setPosition([xPatch(1) yPatch(1) xPatch(2)-xPatch(1) yPatch(3)-yPatch(2)]); % xy
handles.plotWindow.ROIrect(2).setPosition([xPatch(1) zPatch(1) xPatch(2)-xPatch(1) zPatch(3)-zPatch(2)]); % xz
handles.plotWindow.ROIrect(3).setPosition([yPatch(1) zPatch(1) yPatch(3)-yPatch(2) zPatch(3)-zPatch(2)]); % yz
% update figure limits with xy-information
set(handles.plotWindow.hxy,'xlim',(xPatch(2)-xPatch(1))*[-0.25 1.25],'ylim',(yPatch(3)-yPatch(2))*[-0.25 1.25]);
axes(handles.plotWindow.hxy);axis equal

set(handles.plotWindow.hxz,'xlim',(xPatch(2)-xPatch(1))*[-0.25 1.25]);
axes(handles.plotWindow.hxz);axis equal

xz_zlim=get(handles.plotWindow.hxz,'ylim');
set(handles.plotWindow.hyz,'xlim',(yPatch(3)-yPatch(2))*[-0.25 1.25],'ylim',xz_zlim);
axes(handles.plotWindow.hyz);axis equal

% recompute the visible trajectory data
X=handles.plotWindow.xTrj;
Y=camA*vxl*X+camB*ones(1,size(X,2));

set(handles.plotWindow.dat(1),'xdata',Y(1,:),'ydata',Y(2,:)); % xy
set(handles.plotWindow.dat(2),'xdata',Y(1,:),'ydata',Y(3,:)); % xz
set(handles.plotWindow.dat(3),'xdata',Y(2,:),'ydata',Y(3,:)); % yz
function rotate_ROImid(hObject,xy_xz_yzFlag,pmRot)
% hObject       some object handle that can retrieve guidata
% xy_xz_yzFlag  =1,2,3, to indicate which plane to rotate in
% pmRot         rotation direction (+- 1)

handles=guidata(hObject);
dv=pmRot*get(handles.str_angle_increment,'Value')*pi/180; % rotation angle, in radians

switch xy_xz_yzFlag % which rotation matrix to construct
    case 1
        Rdv=[cos(dv) -sin(dv) 0;
             sin(dv)  cos(dv) 0;
            0           0     1];
    case 2
        Rdv=[cos(dv) 0  -sin(dv);
            0        1      0   ;
            sin(dv)  0   cos(dv)];
    case 3
        Rdv=[1      0         0;
             0 cos(dv) -sin(dv);
             0 sin(dv)  cos(dv)];
    otherwise
        error('rotate_ROImid called with undefined flag (only xy_xz_yzFlag  =1,2,3 is defined')
end

% read current ROI parameters
ph=guidata(handles.parentGUI_hObject);
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
px=Data(1);
Nx=Data(2);
Ny=Data(3);
camB=get(ph.camera_b,'Data');
camA=get(ph.camera_A,'Data');

% midpoint of ROI
Ymid=px/2*[Nx+1;Ny+1;0]; 

% rotation transformation
camA=Rdv*camA;
camB=Ymid+Rdv*(camB-Ymid);
set(ph.camera_b,'Data',camB);
set(ph.camera_A,'Data',camA);

% store handles and update trajectory transformation
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);
updateROI_fromMainGUI(hObject);
%% callback functions for data rotation
function str_angle_increment_Callback(hObject, eventdata, handles)
% hObject    handle to str_angle_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of str_angle_increment as text
%        str2double(get(hObject,'String')) returns contents of str_angle_increment as a double

newInc=str2double(get(hObject,'String'));
oldInc=get(hObject,'Value');
if(isnan(newInc) || newInc <=0 )
   errordlg('Angle increment should be a positive number.')
   set(hObject,'String',num2str(oldInc));
else
   set(hObject,'Value',newInc);
end
function str_angle_increment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to str_angle_increment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% rotation callback functions all call rotate_ROImid
function rotate_xy_plus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,1,1);
function rotate_xy_minus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,1,-1);
function rotate_xz_plus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,2,1);
function rotate_xz_minus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,2,-1);
function rotate_yz_plus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,3,1);
function rotate_yz_minus_Callback(hObject, eventdata, handles)
rotate_ROImid(hObject,3,-1);
%% calback for data display update
function read_trj_data_Callback(hObject, eventdata, handles)
% hObject    handle to read_trj_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateTrjData(hObject)
%% callback functions for ROI translations
function translate_xy_Callback(hObject, eventdata, handles)
% hObject    handle to translate_xy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

% read current rectangle positions
rxy=handles.plotWindow.ROIrect(1).getPosition;
rxz=handles.plotWindow.ROIrect(2).getPosition;
ryz=handles.plotWindow.ROIrect(3).getPosition;

% update B: ROI xy origin is at px/2*[1 1 0], to make the origin pixel
% centered at px*[1 1 0].
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
px=Data(1);
camB=get(ph.camera_b,'Data');
camB(1)=camB(1)-rxy(1)+px/2;
camB(2)=camB(2)-rxy(2)+px/2;
set(ph.camera_b,'Data',camB);


% move rectangles back to the relevant origin
rxy(1:2)=px/2*[1 1];
rxz(1)=px/2;
ryz(1)=px/2;
handles.plotWindow.ROIrect(1).setPosition(rxy);    
handles.plotWindow.ROIrect(2).setPosition(rxz);    
handles.plotWindow.ROIrect(3).setPosition(ryz);

% store handles and update trajectory transformation
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);
updateROI_fromMainGUI(hObject);
function translate_xz_Callback(hObject, eventdata, handles)
% hObject    handle to translate_xz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

% read current rectangle positions
rxy=handles.plotWindow.ROIrect(1).getPosition;
rxz=handles.plotWindow.ROIrect(2).getPosition;
ryz=handles.plotWindow.ROIrect(3).getPosition;

% update B: ROI xy origin is at px/2*[-1 -1 0]
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
px=Data(1);
camB=get(ph.camera_b,'Data')
camB(1)=camB(1)-rxz(1)-px/2;
camB(3)=camB(3)-rxz(2);
set(ph.camera_b,'Data',camB);

% move rectangles back to the relevant origin
rxy(1)=-px/2;
rxz(1:2)=[-px/2 0];
ryz(2)=0;
handles.plotWindow.ROIrect(1).setPosition(rxy);    
handles.plotWindow.ROIrect(2).setPosition(rxz);    
handles.plotWindow.ROIrect(3).setPosition(ryz);

% store handles and update trajectory transformation
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);

updateROI_fromMainGUI(hObject);
function translate_yz_Callback(hObject, eventdata, handles)
% hObject    handle to translate_yz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

% read current rectangle positions
rxy=handles.plotWindow.ROIrect(1).getPosition;
rxz=handles.plotWindow.ROIrect(2).getPosition;
ryz=handles.plotWindow.ROIrect(3).getPosition;

% update B: ROI xy origin is at px/2*[-1 -1 0]
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
px=Data(1);
camB=get(ph.camera_b,'Data');
camB(2)=camB(2)-ryz(1)-px/2;
camB(3)=camB(3)-ryz(2);
set(ph.camera_b,'Data',camB);

% move rectangles back to the relevant origin
rxy(2)=-px/2;
rxz(2)=0;
ryz(1:2)=[-px/2 0];
handles.plotWindow.ROIrect(1).setPosition(rxy);    
handles.plotWindow.ROIrect(2).setPosition(rxz);    
handles.plotWindow.ROIrect(3).setPosition(ryz);

% store handles and update trajectory transformation
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);
updateROI_fromMainGUI(hObject);
function translate_data_cm_Callback(hObject, eventdata, handles)
% hObject    handle to translate_data_cm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

% read current rectangle positions
vxl=str2double(get(ph.trj_voxelSize,'String'));
Data=get(ph.camera_pixLength_xRange_yRange,'Data');
px=Data(1);
Nx=Data(2);
Ny=Data(3);
camB=get(ph.camera_b,'Data');
camA=get(ph.camera_A,'Data');

% midpoint of ROI
Ymid=px/2*[Nx+1;Ny+1;0]; 

% c.m. of visible trajectory data
X=handles.plotWindow.xTrj;
Ymean=mean(camA*vxl*X+camB*ones(1,size(X,2)),2); % data c.m. in ROI frame

% compute translation
camB=camB+Ymid-Ymean;
set(ph.camera_b,'Data',camB);

% store handles and update trajectory transformation
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);
updateROI_fromMainGUI(hObject);
%% close and reset callbacks
function reset_A_b_Callback(hObject, eventdata, handles)
% hObject    handle to reset_A_b (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(hObject);
ph=guidata(handles.parentGUI_hObject);

set(ph.camera_b,'Data',zeros(3,1));
set(ph.camera_A,'Data',eye(3));
guidata(hObject,handles);
guidata(handles.parentGUI_hObject,ph);

updateROI_fromMainGUI(hObject);
function close_guiAB_Callback(hObject, eventdata, handles)
% hObject    handle to close_guiAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
palantirAB_gui_CloseFcn([],[],hObject)
