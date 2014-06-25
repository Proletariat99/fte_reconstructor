function varargout = dveGUI_load(varargin)
% DVEGUI_LOAD MATLAB code for dveGUI_load.fig
%      DVEGUI_LOAD, by itself, creates a new DVEGUI_LOAD or raises the existing
%      singleton*.
%
%      H = DVE_LOAD_DATA returns the handle to a new DVE_LOAD_DATA or the handle to
%      the existing singleton*.
%
%      DVE_LOAD_DATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DVE_LOAD_DATA.M with the given input arguments.
%
%      DVE_LOAD_DATA('Property','Value',...) creates a new DVE_LOAD_DATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dveGUI_load_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dveGUI_load_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dveGUI_load

% Last Modified by GUIDE v2.5 03-Aug-2011 12:15:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dveGUI_load_OpeningFcn, ...
  'gui_OutputFcn',  @dveGUI_load_OutputFcn, ...
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


% --- Executes just before dveGUI_load is made visible.
function dveGUI_load_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dveGUI_load (see VARARGIN)

% Choose default command line output for dveGUI_load
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dveGUI_load wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dveGUI_load_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in pm_doy.
function pm_doy_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.doy = contents{get(hObject,'Value')};
disp(['handles.doy from callbackFcn= ', handles.doy])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pm_doy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.doy = contents{get(hObject,'Value')};
disp(['handles.doy from createFcn= ', handles.doy])


function etxt_ti_Callback(hObject, eventdata, handles)
% hObject    handle to etxt_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etxt_ti as text
%        str2double(get(hObject,'String')) returns contents of etxt_ti as a double
handles.ti=get(hObject,'String');
disp(['handles.ti from callbackFcn= ',handles.ti])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function etxt_ti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etxt_ti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.ti=get(hObject,'String');
disp(['handles.ti from createFcn= ',handles.ti])


function etxt_tf_Callback(hObject, eventdata, handles)
handles.tf=get(hObject,'String');
disp(['handles.tf from callbackFcn= ',handles.tf])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function etxt_tf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.tf=get(hObject,'String');
disp(['handles.tf from createFcn= ',handles.tf])


% --- Executes on selection change in pm_sc.
function pm_sc_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.sc = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.sc from callbackFcn= ',handles.sc])


% --- Executes during object creation, after setting all properties.
function pm_sc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.sc = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.sc from createFcn= ',handles.sc])


% --- Executes on selection change in pm_src.
function pm_src_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.src = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.src from callbackFcn= ',handles.src])


% --- Executes during object creation, after setting all properties.
function pm_src_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.src = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.src from createFcn= ',handles.src])


function etxt_nev_Callback(hObject, eventdata, handles)
handles.nev=str2double(get(hObject,'String'));
disp(['handles.nev from callbackFcn= ',num2str(handles.nev)])
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function etxt_nev_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.nev=str2double(get(hObject,'String'));
disp(['handles.nev from CreateFcn= ',num2str(handles.nev)])
guidata(hObject, handles);


function etxt_w_Callback(hObject, eventdata, handles)
handles.w=str2double(get(hObject,'String'));
disp(['handles.w from callbackFcn= ',num2str(handles.w)])
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function etxt_w_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.w=str2double(get(hObject,'String'));
disp(['handles.w from CreateFcn= ',num2str(handles.w)])
guidata(hObject, handles);

% --- Executes on button press in btn_load.
function btn_load_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
handles.dhead=1;
handles.tcf = 24*3600.;  %time conversion factor
handles.yr = sscanf(handles.doy,'%4c%*6c',2);
handles.month = sscanf(handles.doy, '%*5c%2c%*3c',2);
handles.day = sscanf(handles.doy, '%*8c%2c',2);
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);
% File 1 %
fn = ['th',handles.sc,'_',num2str(strcat(handles.yr,handles.month,handles.day)),'_',num2str(strcat(HHi,MMi)),'_',num2str(strcat(HHf,MMf)),'_GSM_v4_',handles.src,'.dat'];
handles.a = importdata(fn,' ',handles.dhead);
disp(['handles.a is loaded.  Class(handles.a.data)  = ',class(handles.a.data)])
disp(['                      Size(handles.a.data)  = ', num2str(size(handles.a.data))])

% File 2 %
handles.n=113;
fnsh=['/FTE/2011/findings/themis_fte_iaga_2011_msheath.txt'];
mshfid = fopen(fnsh);
handles.msh = textscan(mshfid, '%s %s %s %s %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', handles.n);
fclose(mshfid);
disp(['handles.msh is loaded.  Class(handles.msh.data)  = ',class(handles.msh)])
disp(['                      Size(handles.msh.data)  = ', num2str(size(handles.msh))])

% File 3 %
ff='/FTE/2011/findings/themis_fte_iaga_2011_clean_final_timeorder_recons.txt';
fid = fopen(ff);
handles.f = textscan(fid, '%2c %s %s %s %f %s %s %f %f %f', handles.n);
fclose(fid);
disp(['handles.f is loaded.  Class(handles.f.data)  = ',class(handles.f)])
disp(['                      Size(handles.f.data)  = ', num2str(size(handles.f))])

handles.nden=1e6;         %factor for density
handles.nt  =1e6;         %factor for temperature
handles.nv  =1e3;         %factor for velocity
handles.nb  =1e-9;        %factor for magnetic field
handles.kb  =1.38*1e-23;  %Boltzmann's constant
handles.kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
handles.miu =4.0*pi*1e-7; %permeability
handles.mp  =1.67*1e-27;  %proton mass

msh_vx = handles.msh{13}(handles.nev);
msh_vy = handles.msh{14}(handles.nev);
msh_vz = handles.msh{15}(handles.nev);
msh_v = [msh_vx,msh_vy,msh_vz];

disp(['msh_vx=',num2str(msh_vx)]);
disp(['msh_vy=',num2str(msh_vy)]);
disp(['msh_vz=',num2str(msh_vz)]);

handles.V0=msh_v;
nV0=handles.V0/norm(handles.V0);
handles.xs = -nV0;

disp(['h.xs= ',num2str(handles.xs)]);
disp(['h.V0= ',num2str(handles.V0)]);

BX=handles.a.data(:,2); 
BY=handles.a.data(:,3); 
BZ=handles.a.data(:,4); 	% Assigns Bx, By, Bz (from file - GSM)
VX=handles.a.data(:,6); 
VY=handles.a.data(:,7); 
VZ=handles.a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
Nx=handles.a.data(:,10); 					% assigns plasma N density from file (cm^-3)
Tx=handles.a.data(:,11); 					% assigns temp from file (K)
idata=1:length(BX);                 % index from 1 to (# of BX)
%============%
% TIME STUFF %
%============%
sdoy = strcat(num2str(handles.yr),'-',num2str(handles.month),'-',num2str(handles.day));
tstart = strcat(sdoy,'/',handles.ti);
tstop = strcat(sdoy,'/',handles.tf);
%-------------------------------------------------%
% i & f consider only the start/stop time of the  %
% entire data set (on the scale of hours)         %
%-------------------------------------------------%
tsi = datestr(tstart, 'mmmm/dd/yyyy HH:MM:SS:FFF');
tsf = datestr(tstop, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vtsi = datevec(tsi, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vtsf = datevec(tsf, 'mmmm/dd/yyyy HH:MM:SS:FFF');
dni = datenum(vtsi);
dnf = datenum(vtsf);
%---------------------------------%
% 1 & 2 consider stop/start times %
% for FTE from the file or array  %
%---------------------------------%
t1 = handles.f{2}(handles.nev);
t2 = handles.f{4}(handles.nev);
fvts1 = datevec(t1, 'yyyy-mm-dd/HH:MM:SS');
fvts2 = datevec(t2, 'yyyy-mm-dd/HH:MM:SS');
fdn1 = datenum(fvts1);
fdn2 = datenum(fvts2);
%------------------------------------%
% a & b consider adjusted FTE window %
% set by the user input 'w'          %
%------------------------------------%
t1f1 = etime(fvts1,vtsi);   % returns number of seconds from start of window to start of FTE
t1f2 = etime(fvts2,vtsi);   % elapsed time from start of window to stop of FTE (s)
tstep = handles.a.data(:,1);
j = tstep >= t1f1 & tstep < t1f2;
k = find(j);
ns1 = k(1)-handles.w;                % adds two data points to either side
ns2 = k(end) + handles.w;            % of FTE
i1 = ns1;
i2 = ns2;
%================%
% TIME STUFF END %
%================%
handles.bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];				% re-assigns B for some reason, should be the same as above.
handles.vp=[VX(i1:i2) VY(i1:i2) VZ(i1:i2)];
handles.Nx = Nx;
handles.Tx = Tx;
handles.lenb = length(handles.bc(:,1));
lenb = handles.lenb;
% Density for window ta to tb
[EHT,EC,HTR,HTslope,vht]=HTcoef(handles.bc,handles.vp,lenb);
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(handles.bc,lenb);
xtemp = -Kbint;
ytemp = Kbmax;
ztemp = Kbmin;

nxorig=ztemp(1);
nyorig=ztemp(2);
nzorig=ztemp(3);

handles.mvab_theta=acosd(nzorig);
handles.mvab_phi=atand(nyorig/nxorig);  
  
if (nxorig < 0.0) && (nyorig >= 0.0);
    handles.mvab_phi=-handles.mvab_phi+pi;
end
if (nxorig < 0.0) && (nyorig < 0.0); 
    handles.mvab_phi= handles.mvab_phi+pi;
end
if (nxorig >= 0.0) && (nyorig < 0.0); 
    handles.mvab_phi=-handles.mvab_phi+pi*2.;
end

disp(['handles.mvab_phi=',num2str(handles.mvab_phi)])
disp(['handles.mvab_theta=',num2str(handles.mvab_theta)])


handles.mvab_theta=acosd(ztemp);
handles.mvab_phi=atand(ytemp/xtemp);  


dveGUI_controller(handles)
