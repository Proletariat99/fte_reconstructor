function varargout = dve_controller(varargin)
% DVE_CONTROLLER MATLAB code for dve_controller.fig
%      DVE_CONTROLLER, by itself, creates a new DVE_CONTROLLER or raises the existing
%      singleton*.
%
%      H = DVE_CONTROLLER returns the handle to a new DVE_CONTROLLER or the handle to
%      the existing singleton*.
%
%      DVE_CONTROLLER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DVE_CONTROLLER.M with the given input arguments.
%
%      DVE_CONTROLLER('Property','Value',...) creates a new DVE_CONTROLLER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dve_controller_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dve_controller_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_controller

% Last Modified by GUIDE v2.5 02-Aug-2011 11:52:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dve_controller_OpeningFcn, ...
  'gui_OutputFcn',  @dve_controller_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT


% --- Executes just before dve_controller is made visible.
function dve_controller_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dve_controller (see VARARGIN)
handles.output = hObject;  % stolen from controlsuite example
n=113;
fnsh=['/FTE/2011/findings/themis_fte_iaga_2011_msheath.txt'];
mshfid = fopen(fnsh);
msh = textscan(mshfid, '%s %s %s %s %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', n);
fclose(mshfid);

ff='/FTE/2011/findings/themis_fte_iaga_2011_clean_final_timeorder_recons.txt';
fid = fopen(ff);
f = textscan(fid, '%2c %s %s %s %f %s %s %f %f %f', n);
fclose(fid);

handles.msh=msh;
handles.f=f;
handles.nden=1e6;         %factor for density
handles.nt  =1e6;         %factor for temperature
handles.nv  =1e3;         %factor for velocity
handles.nb  =1e-9;        %factor for magnetic field
handles.kb  =1.38*1e-23;  %Boltzmann's constant
handles.kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
handles.miu =4.0*pi*1e-7; %permeability
handles.mp  =1.67*1e-27;  %proton mass
handles.dhead = 1;        % # of header lines in time_step file

% Choose default command line output for dve_controller
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dve_controller wait for user response (see UIRESUME)
% uiwait(handles.Controller);
end

% --- Outputs from this function are returned to the command line.
function varargout = dve_controller_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in btn_thetaup.
function btn_thetaup_Callback(hObject, eventdata, handles)
handles.theta0 = handles.theta0 + handles.pol;
guidata(hObject, handles);
%   function pop_theta(handles.theta0)
%   end
end

function etxt_phi_Callback(hObject, eventdata, handles)

%   function pop_phi(handles.phi0)
%   end

guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
end

% --- Executes during object creation, after setting all properties.
function etxt_phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

function etxt_theta_Callback(hObject, eventdata, handles)
handles.theta0=str2double(get(hObject,'String'));
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

function etxt_azi_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_azi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

function etxt_pol_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_pol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

function btn_phidn_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

function btn_phiup_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end
% --- Executes on button press in btn_thetadn.
function btn_thetadn_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

% --- Executes on button press in btn_display.
function btn_display_Callback(hObject, eventdata, handles)
msh_vx = handles.msh{13}(handles.nev);
msh_vy = handles.msh{14}(handles.nev);
msh_vz = handles.msh{15}(handles.nev);
msh_v = [msh_vx,msh_vy,msh_vz];
V0=msh_v;
nV0=V0/norm(V0);
xs = -nV0;
disp(['xs (from nev) = ', num2str(xs)]);
guidata(hObject, handles);

dhead = handles.dhead;
sc = handles.sc;
src = handles.src;
day = handles.day;
month = handles.month;
yr = handles.yr;
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);

secs = handles.a.data(:,1);
BX=a.data(:,2); BY=a.data(:,3); BZ=a.data(:,4); 	% Assigns Bx, By, Bz (from file - GSM)
VX=a.data(:,6); VY=a.data(:,7); VZ=a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
Nx=a.data(:,10); 					% assigns plasma N density from file (cm^-3)
Tx=a.data(:,11); 					% assigns temp from file (K)
idata=1:length(BX);
guidata(hObject, handles);

disp(['doy=', handles.doy]);
disp(['ti=', handles.ti]);
disp(['tf=', handles.tf]);
disp(['sc=', handles.sc]);
disp(['src=', handles.src]);
disp(['w=', handles.w]);
disp(['nev=', handles.nev]);
disp(['theta0=', handles.theta0]);
disp(['phi0=', handles.phi0]);
disp(['pol=', handles.pol]);
disp(['azi=', handles.azi]);
end

% --- Executes on button press in btn_quit.
function btn_quit_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

% --- Executes on button press in btn_back.
function btn_back_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
end

% --- Executes on button press in btn_load.
function btn_load_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
sc = handles.sc;
src = handles.src;
dhead = handles.dhead;
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);
fn = ['th',sc,'_',num2str(strcat(handles.yr,handles.month,handles.day)),'_',num2str(strcat(HHi,MMi)),'_',num2str(strcat(HHf,MMf)),'_GSM_v4_',src,'.dat'];
handles.a = importdata(fn,' ',dhead);
disp(['handles.a imported successfully' ])
end

function etxt_ti_Callback(hObject, eventdata, handles)
handles.ti=get(hObject,'String');
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
disp(['HHi= ', num2str(HHi)]);
disp(['MMi= ', num2str(MMi)]);
disp(['SSi= ', num2str(SSi)]);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_ti_CreateFcn(hObject, eventdata, handles)
a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end


function etxt_tf_Callback(hObject, eventdata, handles)
handles.tf=get(hObject,'String');
% tf=handles.tf
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);
disp(['HHf= ', num2str(HHf)]);
disp(['MMf= ', num2str(MMf)]);
disp(['SSf= ', num2str(SSf)]);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_tf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in pm_sc.
function pm_sc_Callback(hObject, eventdata, handles)
contents=cellstr(get(hObject,'String'));
handles.sc=contents{get(hObject,'Value')};
disp(['sc= ',handles.sc]);
sc = handles.sc;
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function pm_sc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end


function etxt_nev_Callback(hObject, eventdata, handles)
handles.nev = str2double(get(hObject,'String'));
% disp(['handles.nev= ',num2str(handles.nev)]); % for troubleshooting
end


function etxt_nev_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end


function pm_src_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.src=contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['src= ',handles.src])
end

function pm_src_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end


function etxt_w_Callback(hObject, eventdata, handles)
% hObject    handle to etxt_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etxt_w as text
%        str2double(get(hObject,'String')) returns contents of etxt_w as a double
handles.w=str2double(get(hObject,'String'));
disp(['w= ', handles.w]);
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function etxt_w_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etxt_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

% --- Executes on selection change in pm_doy.
function pm_doy_Callback(hObject, eventdata, handles)
% hObject    handle to pm_doy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pm_doy contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pm_doy
contents = cellstr(get(hObject,'String'));
handles.doy = contents{get(hObject,'Value')};
% disp(handles.doy);
doy = handles.doy;
handles.tcf = 24*3600.;  %time conversion factor
handles.yr = sscanf(doy,'%4c%*6c',2);
handles.month = sscanf(doy, '%*5c%2c%*3c',2);
handles.day = sscanf(doy, '%*8c%2c',2);
disp(num2str([handles.yr, handles.month, handles.day]))
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function pm_doy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pm_doy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
end

% function txt_phi_CreateFcn(phiObject, eventdata, handles)
%   set(phiObject,'String',handles.phi0);
% end

end
end
end
