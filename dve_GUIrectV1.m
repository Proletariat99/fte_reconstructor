function varargout = dve_GUIrectV1(varargin)
%DVE_GUISPIN_V1 M-file for dve_GUIspin_v1.fig
%      DVE_GUISPIN_V1, by itself, creates a new DVE_GUISPIN_V1 or raises the existing
%      singleton*.
%
%      H = DVE_GUISPIN_V1 returns the handle to a new DVE_GUISPIN_V1 or the handle to the
%      existing singleton*.
%
%      DVE_GUISPIN_V1('Property','Value',...) creates a new DVE_GUISPIN_V1 using the given
%      property value pairs. Unrecognized properties are passed via varargin to
%      dve_GUIspin_v1_OpeningFcn.  This calling syntax produces a warning when there is an
%      existing singleton*.
%
%      DVE_GUISPIN_V1('CALLBACK') and DVE_GUISPIN_V1('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVE_GUISPIN_V1.M with the given input arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one instance to
%      run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_GUIspin_v1

% Last Modified by GUIDE v2.5 15-Aug-2011 17:00:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dve_GUIrectV1_OpeningFcn, ...
  'gui_OutputFcn',  @dve_GUIrectV1_OutputFcn, ...
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

function dve_GUIrectV1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn. hObject    handle to figure eventdata
% reserved - to be defined in a future version of MATLAB handles    structure with handles
% and user data (see GUIDATA) varargin   unrecognized PropertyName/PropertyValue pairs
% from the
%            command line (see VARARGIN)
% clear handles; handles = guihandles; doycontents =
% cellstr(get(handles.pm_doy,'String')); handles.doy =
% doycontents{get(handles.pm_doy,'Value')}; handles.ti = get(handles.etxt_ti,'String');
% handles.tf = get(handles.etxt_tf,'String'); sccontents =
% cellstr(get(handles.pm_sc,'String')); handles.sc =
% sccontents{get(handles.pm_sc,'Value')}; srccontents =
% cellstr(get(handles.pm_src,'String')); handles.src =
% srccontents{get(handles.pm_src,'Value')}; handles.w = get(handles.etxt_w,'String');
% handles.nev = get(handles.etxt_nev,'String'); handles.ny =
% get(handles.etxt_ny,'String');
handles.nden=1e6;         %factor for density
handles.nt  =1e6;         %factor for temperature
handles.nv  =1e3;         %factor for velocity
handles.nb  =1e-9;        %factor for magnetic field
handles.kb  =1.38*1e-23;  %Boltzmann's constant
handles.kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
handles.miu =4.0*pi*1e-7; %permeability
handles.mp  =1.67*1e-27;  %proton mass
handles.dhead = 1;
% handles = catstruct(handles,guihandles); disp(['Class of handles= ',class(handles)]);
% disp(['Size of handles= ',num2str(size(handles))]); Choose default command line output
% for dve_GUIspin_v1
handles.output = hObject;
guidata(hObject, handles);

function varargout = dve_GUIrectV1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT); hObject    handle to
% figure eventdata  reserved - to be defined in a future version of MATLAB handles
% structure with handles and user data (see GUIDATA) Get default command line output from
% handles structure
varargout{1} = handles.output;

%--------------------------- Begin GUI Objects code ----------------------------% --
% Initial Data Panel
%--
function pm_doy_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.doy = contents{get(hObject,'Value')};
disp(['handles.doy from callbackFcn= ', handles.doy])
guidata(hObject, handles);
function pm_doy_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.doy = contents{get(hObject,'Value')};
disp(['handles.doy from createFcn= ', handles.doy])
guidata(hObject, handles);
function etxt_ti_Callback(hObject, eventdata, handles)
handles.ti=get(hObject,'String');
disp(['handles.ti from callbackFcn= ',handles.ti])
guidata(hObject, handles);
function etxt_ti_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.ti=get(hObject,'String');
disp(['handles.ti from createFcn= ',handles.ti])
guidata(hObject, handles);
function etxt_tf_Callback(hObject, eventdata, handles)
handles.tf=get(hObject,'String');
disp(['handles.tf from callbackFcn= ',handles.tf])
guidata(hObject, handles);
function etxt_tf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.tf=get(hObject,'String');
disp(['handles.tf from createFcn= ',handles.tf])
guidata(hObject, handles);
function pm_sc_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.sc = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.sc from callbackFcn= ',handles.sc])
function pm_sc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.sc = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.sc from createFcn= ',handles.sc])
function pm_src_Callback(hObject, eventdata, handles)
contents = cellstr(get(hObject,'String'));
handles.src = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.src from callbackFcn= ',handles.src])
function pm_src_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
contents = cellstr(get(hObject,'String'));
handles.src = contents{get(hObject,'Value')};
guidata(hObject, handles);
disp(['handles.src from createFcn= ',handles.src])
function etxt_w_Callback(hObject, eventdata, handles)
handles.w=str2double(get(hObject,'String'));
disp(['handles.w from callbackFcn= ',num2str(handles.w)])
guidata(hObject, handles);
function etxt_w_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.w=str2double(get(hObject,'String'));
disp(['handles.w from CreateFcn= ',num2str(handles.w)])
guidata(hObject, handles);
function etxt_nev_Callback(hObject, eventdata, handles)
handles.nev=str2double(get(hObject,'String'));
disp(['handles.nev from callbackFcn= ',num2str(handles.nev)])
guidata(hObject, handles);
function etxt_nev_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.nev=str2double(get(hObject,'String'));
disp(['handles.nev from CreateFcn= ',num2str(handles.nev)])
guidata(hObject, handles);
%--------------------- Load Data Button Callback -----------------------%
function btn_load_Callback(hObject, eventdata, handles)
handles.dhead=1;
handles.tcf = 24*3600.;  %time conversion factor
handles.nden=1e6;         %factor for density
handles.nt  =1e6;         %factor for temperature
handles.nv  =1e3;         %factor for velocity
handles.nb  =1e-9;        %factor for magnetic field
handles.kb  =1.38*1e-23;  %Boltzmann's constant
handles.kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
handles.miu =4.0*pi*1e-7; %permeability
handles.mp  =1.67*1e-27;  %proton mass
handles.yr = sscanf(handles.doy,'%4c%*6c',2);
handles.month = sscanf(handles.doy, '%*5c%2c%*3c',2);
handles.day = sscanf(handles.doy, '%*8c%2c',2);
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);
disp('===========================Begin LOAD Button Data=================================')
guidata(hObject, handles);
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
handles.secs = handles.a.data(:,1);
handles.BX=handles.a.data(:,2);
handles.BY=handles.a.data(:,3);
handles.BZ=handles.a.data(:,4); 	% Assigns handles.BX, By, Bz (from file - GSM)
msh_vx = handles.msh{13}(handles.nev);
msh_vy = handles.msh{14}(handles.nev);
msh_vz = handles.msh{15}(handles.nev);
msh_v = [msh_vx,msh_vy,msh_vz];
handles.V0=msh_v;
handles.nV0=handles.V0/norm(handles.V0);
handles.xsv0 = -handles.nV0;
VX=handles.a.data(:,6);
VY=handles.a.data(:,7);
VZ=handles.a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
handles.Nx=handles.a.data(:,10); 					% assigns plasma N density from file (cm^-3)
handles.Tx=handles.a.data(:,11); 					% assigns temp from file (K)
idata=1:length(handles.BX);
%============%
% TIME STUFF %
%============%
sdoy = strcat(num2str(handles.yr),'-',num2str(handles.month),'-',num2str(handles.day));
tstart = strcat(sdoy,'/',handles.ti);
tstop = strcat(sdoy,'/',handles.tf);
%-------------------------------------------------%
% i & f consider only the start/stop time of the  % entire data set (on the scale of
% hours)         %
%-------------------------------------------------%
tsi = datestr(tstart, 'mmmm/dd/yyyy HH:MM:SS:FFF');
tsf = datestr(tstop, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vtsi = datevec(tsi, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vtsf = datevec(tsf, 'mmmm/dd/yyyy HH:MM:SS:FFF');
dni = datenum(vtsi);
dnf = datenum(vtsf);
%---------------------------------%
% 1 & 2 consider stop/start times % for FTE from the file or array  %
%---------------------------------%
t1 = handles.f{2}(handles.nev);
t2 = handles.f{4}(handles.nev);
fvts1 = datevec(t1, 'yyyy-mm-dd/HH:MM:SS');
fvts2 = datevec(t2, 'yyyy-mm-dd/HH:MM:SS');
fdn1 = datenum(fvts1);
fdn2 = datenum(fvts2);
%------------------------------------%
% a & b consider adjusted FTE window % set by the user input 'w'          %
%------------------------------------%
t1f1 = etime(fvts1,vtsi);   % returns number of seconds from start of window to start of FTE
t1f2 = etime(fvts2,vtsi);   % elapsed time from start of window to stop of FTE (s)
tstep = handles.a.data(:,1);
j = tstep >= t1f1 & tstep < t1f2;
k = find(j);
ns1 = k(1)-handles.w;                % adds two data points to either side
ns2 = k(end) + handles.w;            % of FTE
handles.i1 = ns1;
handles.i2 = ns2;
%================%
% TIME STUFF END %
%================%

handles.bc=[handles.BX(handles.i1:handles.i2) handles.BY(handles.i1:handles.i2) handles.BZ(handles.i1:handles.i2)];				% re-assigns B for some reason, should be the same as above.
handles.bx = handles.bc(:,1);
handles.by = handles.bc(:,2);
handles.bz = handles.bc(:,3);
handles.ndata = numel(handles.bx);
handles.vp=[VX(handles.i1:handles.i2) VY(handles.i1:handles.i2) VZ(handles.i1:handles.i2)];
handles.lenb = length(handles.bc(:,1));     
% Density for window ta to tb
% [EHT,EC,HTR,HTslope,vht]=HTcoef(handles.bc,handles.vp,handles.lenb);
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(handles.bc,handles.lenb);
handles.Kbmax=Kbmax;
handles.Kbint=Kbint;
handles.Kbmin=Kbmin;
% Basis Vector Notes:   x1=Kbmax=x=x'
%                       x2=Kbint=y=z' x3=Kbmin=z=y' (MP Normal)
handles.ztemp = handles.Kbint;      % 
handles.ytemp = handles.Kbmax;      % 
handles.xtemp = handles.Kbmin;      % MP Normal
handles.nxorig=handles.ztemp(1);
handles.nyorig=handles.ztemp(2);
handles.nzorig=handles.ztemp(3);
handles.mvab_theta=acosd(handles.nzorig);
handles.mvab_phi=atand(abs(handles.nyorig)/abs(handles.nxorig));
if (handles.nxorig < 0.0) && (handles.nyorig >= 0.0);
  handles.mvab_phi=-handles.mvab_phi+180;
end
if (handles.nxorig < 0.0) && (handles.nyorig < 0.0);
  handles.mvab_phi= handles.mvab_phi+180;
end
if (handles.nxorig >= 0.0) && (handles.nyorig < 0.0);
  handles.mvab_phi=-handles.mvab_phi+360;
end
handles.phi=0.0;
handles.theta=0.0;
handles.zs0=handles.Kbint;
disp(['handles.mvab_phi=',num2str(handles.mvab_phi)]);
disp(['handles.mvab_theta=',num2str(handles.mvab_theta)]);
disp(['handles.zs0 = ',num2str(handles.zs0,' %1.3f')]);
set(handles.txt_phi,'String',num2str(handles.phi, '%3.1f'));
set(handles.etxt_phi,'String',num2str(handles.phi, '%3.1f'));
set(handles.txt_theta,'String',num2str(handles.theta, '%3.1f'));
set(handles.etxt_theta,'String',num2str(handles.theta,'%3.1f'));
set(uicontrol(handles.txt_theta0),'String',num2str(handles.mvab_theta, '%3.1f'));
set(uicontrol(handles.txt_phi0),'String',num2str(handles.mvab_phi, '%3.1f'));
set(uicontrol(handles.txt_zs0x),'String',num2str(handles.zs0(1), '%1.2f'));
set(uicontrol(handles.txt_zs0y),'String',num2str(handles.zs0(2), '%1.2f'));
set(uicontrol(handles.txt_zs0z),'String',num2str(handles.zs0(3), '%1.2f'));
disp(['handles.phi=',num2str(handles.phi)]);
disp(['handles.theta=',num2str(handles.theta)]);
guidata(hObject, handles);
disp(strcat('n*orig values = ',num2str(handles.nxorig, ' %1.3f'),',',num2str(handles.nyorig, ' %1.3f'),',',num2str(handles.nzorig, ' %1.3f')))
disp(strcat('V0 (from msh file) = ', num2str(handles.V0,'%3.2f')))
disp('=========================== End Load Button Callback Data =========================')
disp(' ')
disp(' ')
disp(' ')


%--
% Explicit Angle Settings Panel
%--
function etxt_phi_Callback(hObject, eventdata, handles)
handles.phi=str2double(get(hObject,'String'));
disp(['phi from CallbackFcn= ', num2str(handles.phi)])
guidata(hObject, handles);
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
uicontrol(handles.etxt_theta);      %set focus to the next field.
guidata(hObject, handles);          %update again
function etxt_phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.phi=get(hObject,'String');
disp(['handles.phi from CreateFcn= ', num2str(handles.phi)])
function etxt_theta_Callback(hObject, eventdata, handles)
handles.theta=str2double(get(hObject,'String'));
disp(['theta from CallbackFcn= ', num2str(handles.theta)])
guidata(hObject, handles);          %update handles.
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
uicontrol(handles.etxt_azi);
guidata(hObject, handles);          %update again
function etxt_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.theta=str2double(get(hObject,'String'));
disp(['handles.theta from CreateFcn= ', num2str(handles.theta)])
guidata(hObject, handles);
function etxt_azi_Callback(hObject, eventdata, handles)
handles.azi=str2double(get(hObject,'String'));
disp(['handles.azi from CallbackFcn= ', num2str(handles.azi)])
guidata(hObject, handles);
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
uicontrol(handles.etxt_pol);
guidata(hObject, handles);          %update again
function etxt_azi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.azi=str2double(get(hObject,'String'));
disp(['handles.azi from CreateFcn= ', num2str(handles.azi)])
guidata(hObject, handles);
function etxt_pol_Callback(hObject, eventdata, handles)
handles.pol=str2double(get(hObject,'String'));
disp(['handles.pol from CallbackFcn= ', num2str(handles.pol)])
guidata(hObject, handles);
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
uicontrol(handles.btn_display);
guidata(hObject, handles);          %update again
function etxt_pol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.pol=str2double(get(hObject,'String'));
disp(['handles.pol from CreateFcn= ', num2str(handles.pol)])

guidata(hObject, handles);
function etxt_ny_Callback(hObject, eventdata, handles)
handles.ny=str2double(get(hObject,'String'));
disp(['handles.ny from CallbackFcn= ', num2str(handles.ny)])
guidata(hObject, handles);
function etxt_ny_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.ny=str2double(get(hObject,'String'));
disp(['handles.ny from CreateFcn= ', num2str(handles.ny)])
guidata(hObject, handles);
%--
% Manual Angle Adjust Panel
%--
function txt_theta_CreateFcn(thetaObject, eventdata, handles)
handles.txt_theta=str2double(get(thetaObject,'String'));
disp(['handles.txt_theta from CreateFcn= ',num2str(handles.txt_theta)]);
guidata(thetaObject, handles);
function txt_phi_CreateFcn(hObject, eventdata, handles)
handles.txt_phi=str2double(get(hObject,'String'));
disp(['handles.txt_phi from CreateFcn= ',num2str(handles.txt_phi)]);
guidata(hObject, handles);
function btn_phiup_Callback(hObject, eventdata, handles)
handles.phi=handles.phi+handles.azi;
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
set(handles.etxt_phi,'String',num2str(handles.phi));
guidata(hObject, handles);
function btn_phidn_Callback(hObject, eventdata, handles)
handles.phi=handles.phi-handles.azi;
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
set(handles.etxt_phi,'String',num2str(handles.phi));
guidata(hObject, handles);
function btn_thetaup_Callback(hObject, eventdata, handles)
handles.theta=handles.theta+handles.pol;
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
set(handles.etxt_theta,'String',num2str(handles.theta));
guidata(hObject, handles);
function btn_thetadn_Callback(hObject, eventdata, handles)
handles.theta=handles.theta-handles.pol;
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
set(handles.etxt_theta,'String',num2str(handles.theta));
guidata(hObject, handles);
%--
% Control Panel
%--
function btn_clear_Callback(hObject, eventdata, handles)
clear handles;
handles = guihandles;
doycontents = cellstr(get(handles.pm_doy,'String'));
handles.doy = doycontents{get(handles.pm_doy,'Value')};
handles.ti = get(handles.etxt_ti,'String');
handles.tf = get(handles.etxt_tf,'String');
sccontents = cellstr(get(handles.pm_sc,'String'));
handles.sc = sccontents{get(handles.pm_sc,'Value')};
srccontents = cellstr(get(handles.pm_src,'String'));
handles.src = srccontents{get(handles.pm_src,'Value')};
handles.w = get(handles.etxt_w,'String');
handles.nev = get(handles.etxt_nev,'String');
handles.ny =  get(handles.etxt_ny,'String');


function btn_next_Callback(hObject, eventdata, handles)
function btn_exit_Callback(hObject, eventdata, handles)
function btn_save_Callback(hObject, eventdata, handles)
function btn_break_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
disp('remove break or ctrl-option-r to continue')

%----------- Display Button Callback -------------%

function btn_display_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla(handles.axes1);
disp('==================== 1st Coordinate Transform Data ====================')
disp('Transforming zs from GSM to local, using MVAB Kbint value')
vza=handles.V0*handles.zs0';
vzav=vza*handles.zs0;
vhtsv=handles.V0-vzav;
vhts=sqrt(vhtsv*vhtsv');
vhtsn=vhtsv./vhts;
handles.xs0=-vhtsn;
handles.ys0=cross(handles.zs0,handles.xs0);
handles.ys0=handles.ys0./norm(handles.ys0);
disp(' ')
disp(['Starting xs (from MVAB angles) = ', num2str(handles.xs0, ' %3.3f')])
disp(['Starting ys (from MVAB angles) = ', num2str(handles.ys0, ' %3.3f')])
disp(['Starting zs (from MVAB angles) = ', num2str(handles.zs0, ' %3.3f')])
disp(['Dot Product of xs & ys = ', num2str(handles.xs0*handles.ys0')])
disp(['Dot Product of zs & xs = ', num2str(handles.zs0*handles.xs0')])
disp(['Dot Product of ys & zs = ', num2str(handles.ys0*handles.zs0')])
disp('Getting adjustment angles for next rotation...')
disp('==================== End 1st Coordinate Transform  ====================')
disp(' ')






function plotstuff(hObject, eventdata,handles)
%----------------------------------------------------------------------------%
%=============================== SubBlock 1 =================================%
%----------------------------------------------------------------------------%
[Aup,Aupdy,Aupdx]=gsup(handles.c1,handles.x,handles.y,handles.u,handles.udy,handles.udx,handles.fS1);
[Adn,Adndy,Adndx]=gsdn(handles.c1,handles.x,handles.y,handles.u,handles.udy,handles.udx,handles.fS1);
Aup=Aup+Adn;
Aup(handles.mid,:)=Aup(handles.mid,:)/2;
Aupdy=Aupdy+Adndy;
Aupdy(handles.mid,:)=Aupdy(handles.mid,:)/2;
Aupdx=Aupdx+Adndx;
Aupdx(handles.mid,:)=Aupdx(handles.mid,:)/2;
gg=1:handles.ny;
Bzup=zeros(handles.ny,handles.nx);
for j=1:handles.ny
  for i=1:handles.nx
    Bzup(j,i)=polyval(handles.fZ1,Aup(j,i));
    if Aup(j,i)>handles.amax1
      Bzup(j,i)=polyval(handles.fZ1,handles.amax1);
    elseif Aup(j,i)<handles.amin1
      Bzup(j,i)=polyval(handles.fZ1,handles.amin1);
    end
  end
end
B3 = Bzup(gg,:)*handles.b0*1e9;              % defines 3rd axis
pcolor(handles.x*handles.L0,handles.y(gg)*handles.L0,B3);            % displays 3rd axis (bz) in faceted shading
shading interp                       % changes shading to interpolated
minb3=min(min(B3));
maxb3=max(max(B3));
caxis([minb3 maxb3]);
zbar=colorbar('vertical');
set(zbar,'Fontsize',8);
set(get(zbar,'XLabel'),'String','Bz [nT]','Rotation',0,'HorizontalAlignment','left','Fontsize',10);
hold on;
[cc,hh]=contour(handles.x*handles.L0,handles.y(gg)*handles.L0,Aup(gg,:),[-4:0.11:4],'k');
set(hh,'linewidth',1.0);
quiv1=quiver(handles.xa0,zeros(1,handles.ndata),handles.bxs',handles.bys',0.3,'w');  % plots white arrows accross middle.
set(quiv1,'linewidth',2.0);                          % sets linewidth for quiver arrows
set(handles.axes1,'fontsize',10,'TickDir','out','linewidth',1.0);
axis equal
handles.xdist = handles.x*handles.L0;
% fudge factor for x range.
if handles.xdist <= 0
  axis([min(handles.x*handles.L0), 0, min(handles.y*handles.L0), max(handles.y*handles.L0)])
else
  axis([0, max(handles.x*handles.L0), min(handles.y*handles.L0), max(handles.y*handles.L0)])
end
xlabel('x [km]','fontsize',12)
ylabel('y [km]','fontsize',12)
%--
% Vector Field Panel
%--
function etxt_fieldfactor_Callback(hObject, eventdata, handles)
handles.etxt_ff=str2double(get(hObject,'String'));
disp(['etxt_ff from Callback=',num2str(handles.etxt_ff)])
guidata(hObject,handles);
function etxt_fieldfactor_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.etxt_ff=str2double(get(hObject,'String'));
disp(['etxt_ff from CreateFcn=',num2str(handles.etxt_ff)])
guidata(hObject,handles);

%------------------- Vector Load Button ------------------------------------
function btn_vectorload_Callback(hObject, eventdata, handles)
