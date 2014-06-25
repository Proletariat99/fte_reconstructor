function varargout = dve_GUIfinal(varargin)
%DVE_GUIFINAL M-file for dve_GUIfinal.fig
%      DVE_GUIFINAL, by itself, creates a new DVE_GUIFINAL or raises the existing
%      singleton*.
%
%      H = DVE_GUIFINAL returns the handle to a new DVE_GUIFINAL or the handle to
%      the existing singleton*.
%
%      DVE_GUIFINAL('Property','Value',...) creates a new DVE_GUIFINAL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to dve_GUIfinal_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DVE_GUIFINAL('CALLBACK') and DVE_GUIFINAL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVE_GUIFINAL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_GUIfinal

% Last Modified by GUIDE v2.5 09-Aug-2011 14:44:13

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dve_GUIfinal_OpeningFcn, ...
  'gui_OutputFcn',  @dve_GUIfinal_OutputFcn, ...
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


function dve_GUIfinal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)
handles.nden=1e6;         %factor for density
handles.nt  =1e6;         %factor for temperature
handles.nv  =1e3;         %factor for velocity
handles.nb  =1e-9;        %factor for magnetic field
handles.kb  =1.38*1e-23;  %Boltzmann's constant
handles.kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
handles.miu =4.0*pi*1e-7; %permeability
handles.mp  =1.67*1e-27;  %proton mass
handles.dhead = 1;
handles = catstruct(handles,guihandles);
% disp(['Class of handles= ',class(handles)]);
% disp(['Size of handles= ',num2str(size(handles))]);


% Choose default command line output for dve_GUIfinal
handles.output = hObject;

guidata(hObject, handles);

function varargout = dve_GUIfinal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--------------------------- Begin GUI Objects code ----------------------------%
%--
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
% handles = catstruct(handles,guihandles);
handles.dhead=1;
handles.tcf = 24*3600.;  %time conversion factor
handles.yr = sscanf(handles.doy,'%4c%*6c',2);
handles.month = sscanf(handles.doy, '%*5c%2c%*3c',2);
handles.day = sscanf(handles.doy, '%*8c%2c',2);
%dir
% exist
HHi = sscanf(handles.ti,'%2c%0.2*4c');
MMi = sscanf(handles.ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(handles.ti, '%*6c%2c0', 2);
HHf = sscanf(handles.tf,'%2c%0.2*4c');
MMf = sscanf(handles.tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(handles.tf, '%*6c%2c0', 2);
disp('===============Begin LOAD Button Data==================')
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
handles.V0=msh_v;
handles.nV0=handles.V0/norm(handles.V0);
handles.xsv0 = -handles.nV0;
handles.secs = handles.a.data(:,1);
BX=handles.a.data(:,2);
BY=handles.a.data(:,3);
BZ=handles.a.data(:,4); 	% Assigns Bx, By, Bz (from file - GSM)
VX=handles.a.data(:,6);
VY=handles.a.data(:,7);
VZ=handles.a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
handles.Nx=handles.a.data(:,10); 					% assigns plasma N density from file (cm^-3)
handles.Tx=handles.a.data(:,11); 					% assigns temp from file (K)
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
handles.i1 = ns1;
handles.i2 = ns2;
%================%
% TIME STUFF END %
%================%
handles.bc=[BX(handles.i1:handles.i2) BY(handles.i1:handles.i2) BZ(handles.i1:handles.i2)];				% re-assigns B for some reason, should be the same as above.
handles.bx = handles.bc(:,1);
handles.by = handles.bc(:,2);
handles.bz = handles.bc(:,3);
handles.vp=[VX(handles.i1:handles.i2) VY(handles.i1:handles.i2) VZ(handles.i1:handles.i2)];
handles.Nx = handles.Nx;
handles.Tx = handles.Tx;
handles.lenb = length(handles.bc(:,1));
% Density for window ta to tb
[EHT,EC,HTR,HTslope,vht]=HTcoef(handles.bc,handles.vp,handles.lenb);
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(handles.bc,handles.lenb);
xtemp = -Kbint;
ytemp = Kbmax;
ztemp = Kbmin;
nxorig=ztemp(1);
nyorig=ztemp(2);
nzorig=ztemp(3);
handles.mvab_theta=acosd(nzorig);
handles.mvab_phi=atand(nyorig/nxorig);
handles.prelim_zs = Kbmin;
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
handles.phi=handles.mvab_phi;
set(handles.txt_phi,'String',num2str(handles.mvab_phi));
set(handles.etxt_phi,'String',num2str(handles.mvab_phi));
handles.theta=handles.mvab_theta;
set(handles.txt_theta,'String',num2str(handles.mvab_theta));
set(handles.etxt_theta,'String',num2str(handles.mvab_theta));
disp(['Theta=',num2str(handles.theta)])
disp(['Phi=',num2str(handles.phi)])
disp(['Spin=',num2str(handles.spin)])
disp(['actual msh V0 (V0)=',num2str(handles.V0)])
disp(['normalized V0 (nV0)=',num2str(handles.nV0)])
disp(['xs from V0 (-nV0)=',num2str(handles.xsv0)])
disp(['bx(1st row) = ',num2str(handles.bx(1,:))])
disp(['by(1st row) = ',num2str(handles.by(1,:))])
disp(['bz(1st row) = ',num2str(handles.bz(1,:))])
disp(['bx(last row)= ',num2str(handles.bx(end,:))])
disp(['by(last row)= ',num2str(handles.by(end,:))])
disp(['bz(last row)= ',num2str(handles.bz(end,:))])
disp(['bc(1st row) = ',num2str(handles.bc(1,:))])
disp(['bc(last row)= ',num2str(handles.bc(end,:))])
disp(['vp(1st row) = ',num2str(handles.vp(1,:))])
disp(['vp(last row)= ',num2str(handles.vp(end,:))])
disp('--------------')
disp(['prelim_zs (from MVAB)= ',num2str(handles.prelim_zs)])
disp('--------------')
disp(['===============End LOAD Button Data=================='])
guidata(hObject, handles);

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
uicontrol(handles.etxt_spin);
guidata(hObject, handles);          %update again

function etxt_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.theta=str2double(get(hObject,'String'));
disp(['handles.theta from CreateFcn= ', num2str(handles.theta)])
guidata(hObject, handles);

function etxt_spin_Callback(hObject, eventdata, handles)
handles.spin=str2double(get(hObject,'String'));
disp(['spin from CallbackFcn= ', num2str(handles.spin)])
guidata(hObject, handles);
set(uicontrol(handles.txt_spin),'String',num2str(handles.spin));
uicontrol(handles.etxt_azi);
guidata(hObject, handles);          %update again

function etxt_spin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.spin=str2double(get(hObject,'String'));
disp(['handles.spin from CreateFcn= ', num2str(handles.spin)])
guidata(hObject,handles);

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
uicontrol(handles.etxt_axi);
guidata(hObject, handles);          %update again

function etxt_pol_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.pol=str2double(get(hObject,'String'));
disp(['handles.pol from CreateFcn= ', num2str(handles.pol)])

guidata(hObject, handles);

function etxt_axi_Callback(hObject, eventdata, handles)
handles.axi=str2double(get(hObject,'String'));
disp(['handles.axi from CallbackFcn= ', num2str(handles.axi)])
guidata(hObject, handles);
set(uicontrol(handles.txt_spin),'String',num2str(handles.spin));
uicontrol(handles.btn_display);
guidata(hObject, handles);

function etxt_axi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.axi=str2double(get(hObject,'String'));
disp(['handles.axi from CreateFcn= ', num2str(handles.axi)])
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
function txt_theta_CreateFcn(hObject, eventdata, handles)
handles.txt_theta=str2double(get(hObject,'String'));
disp(['handles.txt_theta from CreateFcn= ',num2str(handles.txt_theta)]);
guidata(hObject, handles);

function txt_phi_CreateFcn(hObject, eventdata, handles)
handles.txt_phi=str2double(get(hObject,'String'));
disp(['handles.txt_phi from CreateFcn= ',num2str(handles.txt_phi)]);
guidata(hObject, handles);

function btn_phiup_Callback(hObject, eventdata, handles)
handles.phi=handles.phi+handles.azi;
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
guidata(hObject, handles);

function btn_phidn_Callback(hObject, eventdata, handles)
handles.phi=handles.phi-handles.azi;
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
guidata(hObject, handles);

function btn_thetaup_Callback(hObject, eventdata, handles)
handles.theta=handles.theta+handles.pol;
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
guidata(hObject, handles);

function btn_thetadn_Callback(hObject, eventdata, handles)
handles.theta=handles.theta-handles.pol;
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
guidata(hObject, handles);

function btn_spindn_Callback(hObject, eventdata, handles)
handles.spin=handles.spin-handles.axi;
set(uicontrol(handles.txt_spin),'String',num2str(handles.spin));
guidata(hObject, handles);

function btn_spinup_Callback(hObject, eventdata, handles)
handles.spin=handles.spin+handles.axi;
set(uicontrol(handles.txt_spin),'String',num2str(handles.spin));
guidata(hObject, handles);

%--
% Control Panel
%--
function btn_prev_Callback(hObject, eventdata, handles)

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
% Importing Handles variables to local variables (clunky, but fast);
guidata(hObject, handles);
handles.ndata = handles.lenb;
handles.dTim1 = handles.secs(handles.i1:handles.i2)-handles.secs(handles.i1);               % my code (dd, 2011)
%Calculate z components based on phi/theta
zx = sind(handles.theta)*cosd(handles.phi);
zy = sind(handles.theta)*sind(handles.phi);
zz = cosd(handles.theta);
handles.zs = [zx zy zz];
%   handles.ys = cross(handles.zs,handles.xsv0);
xx = sind(handles.theta+90)*cosd(handles.phi+0);
xy = sind(handles.theta+90)*sind(handles.phi+0);
xz = cosd(handles.theta+90);
handles.xs = [xx xy xz];
handles.ys = cross(handles.zs,handles.xs);

disp(' ')
disp(' ')
disp('ooooooooooooooooooooo Begin Display Button Variables oooooooooooooooooooooooo')
disp(' ')
disp(['theta (beginning of display)=',num2str(handles.theta)])
disp(['phi (beginning of display)=',num2str(handles.phi)])
disp(['handles.zs (from this theta & phi) =', num2str(handles.zs)])
disp(['handles.xs (from this theta & phi) =', num2str(handles.xs)])
disp(['handles.ys (zs cross xs) =', num2str(handles.ys)])
disp('---- ')


% rotate new coordinate axis by "spin" angle
handles.nszs=handles.zs;

handles.spun_xs = handles.xs*cosd(handles.spin)-handles.ys*sind(handles.spin);
handles.spun_ys = handles.xs*sind(handles.spin)+handles.ys*cosd(handles.spin);
handles.nsys = handles.spun_ys/norm(handles.spun_ys);
handles.nsxs = handles.spun_xs/norm(handles.spun_xs);
%   disp(['norm(handles.zs) after spin=',num2str(norm(handles.zs))])
%   disp(['norm(handles.xs) after spin=',num2str(norm(handles.xs))])
%   disp(['norm(handles.ys) after spin=',num2str(norm(handles.ys))])

% project B into reconstruction coord.
% Note: nsxs = 'normalized spun xs'
handles.bxs=(handles.bc.*handles.nb)*handles.nsxs';
handles.bys=(handles.bc.*handles.nb)*handles.nsys';
handles.bzs=(handles.bc.*handles.nb)*handles.nszs';
%   handles.bxs = bxs;
%   handles.bys = bys;
%   handles.bzs = bzs;
handles.vxs = -handles.V0*handles.nsxs';
handles.xa0=handles.vxs.*handles.dTim1;
handles.iix = handles.i1:handles.i2;
handles.dzn=handles.Nx(handles.iix);                    % density reconstruction by matrix
handles.Tp = handles.Tx(handles.iix)/handles.kb2;
%Project V into reconstruction coord.
disp(['-----------B and V in spun coordinate system---------'])
disp(['spun theta =',num2str(handles.theta)])
disp(['spun phi   =',num2str(handles.phi)])
disp(['spun spin  =',num2str(handles.spin)])
disp(['handles.nszs =', num2str(handles.nszs)])
disp(['handles.nsxs =', num2str(handles.nsxs)])
disp(['handles.nsys =', num2str(handles.nsys)])
disp('---- ')
vxc=zeros(1,handles.ndata);
vyc=zeros(1,handles.ndata);
vzc=zeros(1,handles.ndata);
for i=1:handles.ndata
  vxc(i)=(handles.vp(i,:)-handles.V0)*handles.nsxs';
  vyc(i)=(handles.vp(i,:)-handles.V0)*handles.nsys';
  vzc(i)=(handles.vp(i,:)-handles.V0)*handles.nszs';
end
%Calculate vector potential A along y=0
ndA=zeros(1,handles.lenb);
for i=2:handles.lenb                       %
  ndA(i)=-(handles.bys(i)+handles.bys(i-1))*(handles.dTim1(i)-handles.dTim1(i-1))*handles.vxs*handles.nv*0.5;  % Determines A for different
end
A1=cumsum(ndA);                   % with every additional deltat, we add more area
A1=A1';                           % transpose A1
bb=sqrt(handles.bxs.^2+handles.bys.^2+handles.bzs.^2);    % magnitude of in-plane bfield vector
pp=handles.dzn.*handles.Tp*handles.kb*handles.nden;                     % new pp factor (without nt)
bmax=max(bb);                           % maximum b magnitude
nmax=max(handles.dzn)*handles.nden;                     % density maximum value
%Normalized factors
handles.b0=bmax;                                % max b arrow
p0=bmax*bmax/handles.miu;                       % max pressure value
n0=nmax;                                % max density value from 5 lines above.
vv0=sqrt(handles.b0^2/(handles.miu*n0*handles.mp)*1e-6);
A0=max(abs(A1));
handles.L0=(A0/bmax)*1e-3;
disp(['L0=',num2str(handles.L0)])
handles.pbz=pp./p0+((handles.bzs./handles.b0).^2)/2;
handles.An=A1./A0;
guidata(hObject, handles);
handles.amax1=max(handles.An);
handles.amin1=min(handles.An);
handles.fS1=polyfit(handles.An,handles.pbz,4); %Pt(A)
handles.fZ1=polyfit(handles.An,handles.bzs/handles.b0,5); %Bz(A)
%--interpolation--%
handles.nx=handles.ndata+(handles.ndata-1)*3;               % is this just an approximation for distance?
handles.xi=handles.xa0(1):(handles.xa0(handles.ndata)-handles.xa0(1))/(handles.nx-1):handles.xa0(handles.ndata);
bxi=interp1(handles.xa0,handles.bxs./handles.b0,handles.xi,'spline');
byi=interp1(handles.xa0,handles.bys./handles.b0,handles.xi,'spline');
ht=0;
py=0.1/1;
handles.mid=round(handles.ny/2)+ht;
handles.x=handles.xi./handles.L0;
guidata(hObject,handles);
hx=handles.x(2)-handles.x(1);%uniform grids
hy=py*hx;
handles.y = zeros(1,handles.ny);
for j=1:handles.ny
  handles.y(j)=(j-handles.mid)*hy;
end
%---------------------------%
clear('ndA');
ndA(1)=0;
for i=2:handles.nx
  ndA(i)=-(byi(i)+byi(i-1))*(handles.x(i)-handles.x(i-1))*0.5;
end
A2=cumsum(ndA);
handles.u=zeros(handles.ny,handles.nx);
handles.udy=zeros(handles.ny,handles.nx);
handles.udx=zeros(handles.ny,handles.nx);
handles.u(handles.mid,:)=A2;
handles.udy(handles.mid,:)=bxi;
handles.udx(handles.mid,:)=byi;
handles.c1=[handles.nx handles.ny handles.mid hx hy handles.amax1 handles.amin1];
%--plot limits--
handles.ymin = min(handles.y*handles.L0);
handles.ymax = max(handles.y*handles.L0);
handles.xmin = min(handles.x*handles.L0);
handles.xmax = max(handles.x*handles.L0);
handles.hafx = handles.xmax*0.5;
handles.hafy = handles.ymax*0.5;
guidata(hObject, handles);

%-- call plot function
plotstuff(hObject,eventdata,handles);
%-- refresh/plot vector field
guidata(hObject, handles);

phick=get(handles.cb_showphi,'Value');
thtck=get(handles.cb_showtheta,'Value');
spnck=get(handles.cb_showspin,'Value');
ckarr=[phick thtck spnck];
str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)), num2str(ckarr(3)));
switch str_ck
  case '000'
    cla(handles.axes6);
  case '111'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showphi_Callback(hObject,eventdata,handles);
    cb_showtheta_Callback(hObject,eventdata,handles);
    cb_showspin_Callback(hObject,eventdata,handles);
  case '100'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showphi_Callback(hObject,eventdata,handles);
  case '110'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showphi_Callback(hObject,eventdata,handles);
    cb_showtheta_Callback(hObject,eventdata,handles);
  case '011'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showtheta_Callback(hObject,eventdata,handles);
    cb_showspin_Callback(hObject,eventdata,handles);
  case '001'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showspin_Callback(hObject,eventdata,handles);
  case '101'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showspin_Callback(hObject,eventdata,handles);
    cb_showphi_Callback(hObject,eventdata,handles);
  case '010'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showtheta_Callback(hObject,eventdata,handles);
  otherwise
    disp('Some other issue exists. Oh, and go brush up your linear algebra, Captain Clunky.')
end

disp(' ')
disp('oooooooooooooooooo-- End Display Button Data --oooooooooooooooooo')

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
if handles.xdist <= 0
  axis([min(handles.x*handles.L0), 0, min(handles.y*handles.L0), max(handles.y*handles.L0)])
else
  axis([0, max(handles.x*handles.L0), min(handles.y*handles.L0), max(handles.y*handles.L0)])
end
xlabel('x [km]','fontsize',12)
ylabel('y [km]','fontsize',12)
%--
%Curve Fitting Processing
%--
handles.fS1=polyfit(handles.An,handles.pbz,4); %Pt(A)
handles.fZ1=polyfit(handles.An,handles.bzs/handles.b0,5); %Bz(A)
%   amax1=max(handles.An); amin1=min(handles.An);
handles.ax1=handles.amin1:0.01:handles.amax1;
df1=polyder(handles.fS1); %dPt/dA
yys1=polyval(handles.fS1,handles.ax1);
yyz1=polyval(handles.fZ1,handles.ax1);
guidata(hObject,handles);

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


%------------- check boxes (cb) for vector plots --------------%
%---- phi -----%
function cb_showphi_Callback(hObject, eventdata, handles)
if (get(handles.cb_showphi,'Value') == get(handles.cb_showphi,'Max'))
  % Checkbox is checked-take appropriate action
  axes(handles.axes6);
  hold on;
  handles.xdist = handles.x*handles.L0;
  % limits on frame.
  if handles.xdist <= 0
    axis([min(handles.x*handles.L0), 0, min(handles.y*handles.L0), max(handles.y*handles.L0)])
  else
    axis([0, max(handles.x*handles.L0), min(handles.y*handles.L0), max(handles.y*handles.L0)])
  end
  quiv0=quiver(handles.xa0,0,handles.bxs',handles.bys',0.25,'k','linewidth',2);
  % add phi
  handles.quiv1=quiver(handles.xa0,(handles.hafy*0.5),handles.bxp1',handles.byp1',0.25,'b','linewidth', 2);
  handles.quiv2=quiver(handles.xa0,(handles.hafy*1.0),handles.bxp2',handles.byp2',0.25,'b','linewidth', 2);
  handles.quiv3=quiver(handles.xa0,(handles.hafy*1.5),handles.bxp3',handles.byp3',0.25,'b','linewidth', 2);
  % subtract phi
  handles.quiv4=quiver(handles.xa0,-(handles.hafy*0.5),handles.nbxp1',handles.nbyp1',0.25,'b','linewidth', 2);
  handles.quiv5=quiver(handles.xa0,-(handles.hafy*1.0),handles.nbxp2',handles.nbyp2',0.25,'b','linewidth', 2);
  handles.quiv6=quiver(handles.xa0,-(handles.hafy*1.5),handles.nbxp3',handles.nbyp3',0.25,'b','linewidth', 2);
  handles.txt_minusphi=text(0.75*handles.hafx,-handles.hafy*1.75,'MINUS PHI','fontsize',14,'color','b');
  handles.txt_plusphi=text(0.75*handles.hafx,handles.hafy*1.75,'PLUS PHI','fontsize',14,'color','b');
  guidata(hObject,handles);
else
  phick=get(handles.cb_showphi,'Value');
  thtck=get(handles.cb_showtheta,'Value');
  spnck=get(handles.cb_showspin,'Value');
  ckarr=[phick thtck spnck];
  str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)), num2str(ckarr(3)));
  switch str_ck
    case '000'
      cla(handles.axes6);
    case '111'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '100'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
    case '110'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
    case '011'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '001'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
    case '101'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
      cb_showphi_Callback(hObject,eventdata,handles);
    case '010'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
    otherwise
      disp('Some other issue exists. Oh, and go brush up your linear algebra, Captain Clunky.')
  end
end

%------theta------%
function cb_showtheta_Callback(hObject, eventdata, handles)
if (get(handles.cb_showtheta,'Value') == get(handles.cb_showtheta,'Max'))
  axes(handles.axes6);
  hold on;
  handles.xdist = handles.x*handles.L0;
  % limits on frame.
  if handles.xdist <= 0
    axis([min(handles.x*handles.L0), 0, min(handles.y*handles.L0), max(handles.y*handles.L0)])
  else
    axis([0, max(handles.x*handles.L0), min(handles.y*handles.L0), max(handles.y*handles.L0)])
  end
  % add theta (for later use)
  handles.quiv13=quiver(handles.xa0,0,handles.bxs',handles.bys',0.25,'k','linewidth',2);
  handles.quiv7=quiver(handles.xa0,(handles.hafy*0.5),handles.bxt1',handles.byt1',0.25,'g','linewidth', 2);
  handles.quiv8=quiver(handles.xa0,(handles.hafy*1.0),handles.bxt2',handles.byt2',0.25,'g','linewidth', 2);
  handles.quiv9=quiver(handles.xa0,(handles.hafy*1.5),handles.bxt3',handles.byt3',0.25,'g','linewidth', 2);
  % subtract theta
  handles.quiv10=quiver(handles.xa0,-(handles.hafy*0.5),handles.nbxt1',handles.nbyt1',0.25,'g','linewidth', 2);
  handles.quiv11=quiver(handles.xa0,-(handles.hafy*1.0),handles.nbxt2',handles.nbyt2',0.25,'g','linewidth', 2);
  handles.quiv12=quiver(handles.xa0,-(handles.hafy*1.5),handles.nbxt3',handles.nbyt3',0.25,'g','linewidth', 2);
  % text guide
  handles.txt_minustheta=text(0.25*handles.hafx,-handles.hafy*1.75,'MINUS THETA','fontsize',14,'color','g');
  handles.txt_plustheta=text(0.25*handles.hafx,handles.hafy*1.5,'PLUS THETA','fontsize',14,'color','g');
  %   handles.tchildren=get(handles.axes6,'children');
  guidata(hObject,handles);
else
  % Checkbox is not checked-take appropriate action
  cla(handles.axes6);
  phick=get(handles.cb_showphi,'Value');
  thtck=get(handles.cb_showtheta,'Value');
  spnck=get(handles.cb_showspin,'Value');
  ckarr=[phick thtck spnck];
  str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)), num2str(ckarr(3)));
  switch str_ck
    case '000'
      cla(handles.axes6);
    case '111'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '100'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
    case '110'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
    case '011'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '001'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
    case '101'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
      cb_showphi_Callback(hObject,eventdata,handles);
    case '010'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
    otherwise
      disp('Some other issue exists. Oh, and go brush up your linear algebra, Captain Clunky.')
  end
end

%--------spin-------%
function cb_showspin_Callback(hObject, eventdata, handles)
if (get(handles.cb_showspin,'Value') == get(handles.cb_showspin,'Max'))
  axes(handles.axes6);
  hold on;
  handles.xdist = handles.x*handles.L0;
  % limits on frame.
  if handles.xdist <= 0
    axis([min(handles.x*handles.L0), 0, min(handles.y*handles.L0), max(handles.y*handles.L0)])
  else
    axis([0, max(handles.x*handles.L0), min(handles.y*handles.L0), max(handles.y*handles.L0)])
  end
  % add spin
  handles.quiv14=quiver(handles.xa0,0,handles.bxs',handles.bys',0.25,'k','linewidth',2);
  handles.quiv15=quiver(handles.xa0,(handles.hafy*0.5),handles.bxs1',handles.bys1',0.25,'r','linewidth', 2);
  handles.quiv16=quiver(handles.xa0,(handles.hafy*1.0),handles.bxs2',handles.bys2',0.25,'r','linewidth', 2);
  handles.quiv17=quiver(handles.xa0,(handles.hafy*1.5),handles.bxs3',handles.bys3',0.25,'r','linewidth', 2);
  % subtract spin
  handles.quiv18=quiver(handles.xa0,-(handles.hafy*0.5),handles.nbxs1',handles.nbys1',0.25,'r','linewidth', 2);
  handles.quiv19=quiver(handles.xa0,-(handles.hafy*1.0),handles.nbxs2',handles.nbys2',0.25,'r','linewidth', 2);
  handles.quiv20=quiver(handles.xa0,-(handles.hafy*1.5),handles.nbxs3',handles.nbys3',0.25,'r','linewidth', 2);
  % text guide
  handles.txt_minuspin=text(0.25*handles.hafx,-handles.hafy*1.75,'MINUS SPIN','fontsize',14,'color','r');
  hanldes.txt_pluspin=text(0.25*handles.hafx,handles.hafy*1.75,'PLUS SPIN','fontsize',14,'color','r');
  guidata(hObject,handles);
else
  %   % Checkbox is not checked-take appropriate action
  phick=get(handles.cb_showphi,'Value');
  thtck=get(handles.cb_showtheta,'Value');
  spnck=get(handles.cb_showspin,'Value');
  ckarr=[phick thtck spnck];
  str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)), num2str(ckarr(3)));
  switch str_ck
    case '000'
      cla(handles.axes6);
    case '111'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '100'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
    case '110'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
    case '011'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
      cb_showspin_Callback(hObject,eventdata,handles);
    case '001'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
    case '101'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showspin_Callback(hObject,eventdata,handles);
      cb_showphi_Callback(hObject,eventdata,handles);
    case '010'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showtheta_Callback(hObject,eventdata,handles);
    otherwise
      disp('Some other issue exists. Oh, and go brush up your linear algebra, Captain Clunky.')
  end
end

function btn_vectorload_Callback(hObject, eventdata, handles)
gbg_BX=handles.a.data(:,2);
gbg_BY=handles.a.data(:,3);
gbg_BZ=handles.a.data(:,4);
handles.gbg_bc=[gbg_BX(handles.i1:handles.i2) gbg_BY(handles.i1:handles.i2) gbg_BZ(handles.i1:handles.i2)];
handles.theta10 = handles.theta + 1*handles.etxt_ff;
handles.theta20 = handles.theta + 2*handles.etxt_ff;
handles.theta30 = handles.theta + 3*handles.etxt_ff;
handles.phi10 = handles.phi + 1*handles.etxt_ff;
handles.phi20 = handles.phi + 2*handles.etxt_ff;
handles.phi30 = handles.phi + 3*handles.etxt_ff;
handles.spin10 = handles.spin + 1*handles.etxt_ff;
handles.spin20 = handles.spin + 2*handles.etxt_ff;
handles.spin30 = handles.spin + 3*handles.etxt_ff;
%--- zs components by angle
%--- +10 theta
zx10t = sind(handles.theta10)*cosd(handles.phi);
zy10t = sind(handles.theta10)*sind(handles.phi);
zz10t = cosd(handles.theta10);
%--- +20 theta
zx20t = sind(handles.theta20)*cosd(handles.phi);
zy20t = sind(handles.theta20)*sind(handles.phi);
zz20t = cosd(handles.theta20);
%--- +30 theta
zx30t = sind(handles.theta30)*cosd(handles.phi);
zy30t = sind(handles.theta30)*sind(handles.phi);
zz30t = cosd(handles.theta30);
%--- +10 phi
zx10p = sind(handles.theta)*cosd(handles.phi20);
zy10p = sind(handles.theta)*sind(handles.phi20);
zz10p = cosd(handles.theta);
%--- +20 phi
zx20p = sind(handles.theta)*cosd(handles.phi20);
zy20p = sind(handles.theta)*sind(handles.phi20);
zz20p = cosd(handles.theta);
%--- +30 phi
zx30p = sind(handles.theta)*cosd(handles.phi30);
zy30p = sind(handles.theta)*sind(handles.phi30);
zz30p = cosd(handles.theta);
disp(' ')
disp('--------- positive values --------')
disp(['theta used for calculating zs = ',num2str(handles.theta)]);
disp(['  phi used for calculating zs = ',num2str(handles.phi)]);
disp(['theta used for calculating gbg_zs10 = ',num2str(handles.theta10)]);
disp(['  phi used for calculating gbg_zs10 = ',num2str(handles.phi10)]);
disp(['theta used for calculating gbg_zs20 = ',num2str(handles.theta20)]);
disp(['  phi used for calculating gbg_zs20 = ',num2str(handles.phi20)]);
disp(['theta used for calculating gbg_zs30 = ',num2str(handles.theta30)]);
disp(['  phi used for calculating gbg_zs30 = ',num2str(handles.phi30)]);
disp('--------end positive values------')
disp(' ')
guidata(hObject,handles);
%--- xs components by angle (theta)
xx10t = sind(handles.theta10+90)*cosd(handles.phi+0);
xy10t = sind(handles.theta10+90)*sind(handles.phi+0);
xz10t = cosd(handles.theta10+90);
xx20t = sind(handles.theta20+90)*cosd(handles.phi+0);
xy20t = sind(handles.theta20+90)*sind(handles.phi+0);
xz20t = cosd(handles.theta20+90);
xx30t = sind(handles.theta30+90)*cosd(handles.phi+0);
xy30t = sind(handles.theta30+90)*sind(handles.phi+0);
xz30t = cosd(handles.theta30+90);
%--- xs component by angle (phi)
xx10p = sind(handles.theta+90)*cosd(handles.phi10+0);
xy10p = sind(handles.theta+90)*sind(handles.phi10+0);
xz10p = cosd(handles.theta+90);
xx20p = sind(handles.theta+90)*cosd(handles.phi20+0);
xy20p = sind(handles.theta+90)*sind(handles.phi20+0);
xz20p = cosd(handles.theta+90);
xx30p = sind(handles.theta+90)*cosd(handles.phi30+0);
xy30p = sind(handles.theta+90)*sind(handles.phi30+0);
xz30p = cosd(handles.theta+90);
%--- xs component by angle (spin)
xx10s = handles.xs*cosd(handles.spin10)-handles.ys*sind(handles.spin10);
xy10s = handles.xs*sind(handles.spin10)+handles.ys*cosd(handles.spin10);
xz10s = handles.xs(3);
xx20s = handles.xs*cosd(handles.spin20)-handles.ys*sind(handles.spin20);
xy20s = handles.xs*sind(handles.spin20)+handles.ys*cosd(handles.spin20);
xz20s = handles.xs(3);
xx30s = handles.xs*cosd(handles.spin30)-handles.ys*sind(handles.spin30);
xy30s = handles.xs*sind(handles.spin30)+handles.ys*cosd(handles.spin30);
xz30s = handles.xs(3);
% handles.spun_xs = handles.xs*cosd(handles.spin)-handles.ys*sind(handles.spin);
% handles.spun_ys = handles.xs*sind(handles.spin)+handles.ys*cosd(handles.spin);
% handles.xs = [xx xy xz];
% +10 theta zs (garbage)
gbg_zs10t = [zx10t zy10t zz10t];
gbg_xs10t = [xx10t xy10t xz10t];
gbg_ys10t = cross(gbg_zs10t,gbg_xs10t);
% +20 theta zs (garbage)
gbg_zs20t = [zx20t zy20t zz20t];
gbg_xs20t = [xx20t xy20t xz20t];
gbg_ys20t = cross(gbg_zs20t,gbg_xs20t);
% + 30 theta zs (garbage)
gbg_zs30t = [zx30t zy30t zz30t];
gbg_xs30t = [xx10t xy30t xz30t];
gbg_ys30t = cross(gbg_zs30t,gbg_xs30t);
% +10 phi zs (garbage)
gbg_zs10p = [zx10p zy10p zz10p];
gbg_xs10p = [xx10p xy10p xz10p];
gbg_ys10p = cross(gbg_zs10p,gbg_xs10p);
% +20 phi zs (garbage)
gbg_zs20p = [zx20p zy20p zz20p];
gbg_xs20p = [xx20p xy20p xz20p];
gbg_ys20p = cross(gbg_zs20p,gbg_xs20p);
% + 30 phi zs (garbage)
gbg_zs30p = [zx30p zy30p zz30p];
gbg_xs30p = [xx30p xy30p xz30p];
gbg_ys30p = cross(gbg_zs30p,gbg_xs30p);

% --- FOR REFERENCE --- this is how b*s is calculated
% --- ns*s is normalized, spun, x, y and z axis.
% handles.bxs=(handles.bc.*handles.nb)*handles.nsxs';
% handles.bys=(handles.bc.*handles.nb)*handles.nsys';
% handles.bzs=(handles.bc.*handles.nb)*handles.nszs';
% So I need to spin xx10s' before I calculate theta and phi adjustments.
% then I need to do theta/phi adjustments before I calculate spin adjustments.
% shit.

% x & y vectors from spin
handles.bxs1=(handles.bc.*handles.nb)*xx10s';
handles.bxs2=(handles.bc.*handles.nb)*xx20s';
handles.bxs3=(handles.bc.*handles.nb)*xx30s';
handles.bys1=(handles.bc.*handles.nb)*xy10s';
handles.bys2=(handles.bc.*handles.nb)*xy20s';
handles.bys3=(handles.bc.*handles.nb)*xy30s';
% x vectors
handles.bxt1 = (handles.bc.*handles.nb)*gbg_xs10t';
handles.bxt2 = (handles.bc.*handles.nb)*gbg_xs20t';
handles.bxt3 = (handles.bc.*handles.nb)*gbg_xs30t';
handles.bxp1 = handles.bc.*handles.nb*gbg_xs10p';
handles.bxp2 = handles.bc.*handles.nb*gbg_xs20p';
handles.bxp3 = handles.bc.*handles.nb*gbg_xs30p';
% y vectors
handles.byt1 = handles.bc.*handles.nb*gbg_ys10t';
handles.byt2 = handles.bc.*handles.nb*gbg_ys20t';
handles.byt3 = handles.bc.*handles.nb*gbg_ys30t';
handles.byp1 = handles.bc.*handles.nb*gbg_ys10p';
handles.byp2 = handles.bc.*handles.nb*gbg_ys20p';
handles.byp3 = handles.bc.*handles.nb*gbg_ys30p';

disp('-----------first values---------------------------')
disp('---pos x theta---')
disp([' 0deg:',num2str(handles.bxs(1))])
disp(['10deg:',num2str(handles.bxt1(1))])
disp(['20deg:',num2str(handles.bxt2(1))])
disp(['30deg:',num2str(handles.bxt3(1))])
disp('---pos x phi---')
disp([' 0deg:',num2str(handles.bxs(1))])
disp(['10deg:',num2str(handles.bxp1(1))])
disp(['20deg:',num2str(handles.bxp2(1))])
disp(['30deg:',num2str(handles.bxp3(1))])
disp('---pos y theta---')
disp([' 0deg:',num2str(handles.bys(1))])
disp(['10deg:',num2str(handles.byt1(1))])
disp(['20deg:',num2str(handles.byt2(1))])
disp(['30deg:',num2str(handles.byt3(1))])
disp('---pos y phi---')
disp([' 0deg:',num2str(handles.bys(1))])
disp(['10deg:',num2str(handles.byp1(1))])
disp(['20deg:',num2str(handles.byp2(1))])
disp(['30deg:',num2str(handles.byp3(1))])
disp('--------------------------------------------------')
%-----------------NEGATIVE ANGLE ADJUSTMENTS---------%
guidata(hObject,handles);                            %
% negative angles
%--- theta
ntheta10 = handles.theta - 1*handles.etxt_ff;        %
ntheta20 = handles.theta - 2*handles.etxt_ff;        %
ntheta30 = handles.theta - 3*handles.etxt_ff;        %
%--- phi
nphi10 = handles.phi - 1*handles.etxt_ff;            %
nphi20 = handles.phi - 2*handles.etxt_ff;            %
nphi30 = handles.phi - 3*handles.etxt_ff;            %
%--- spin
nspin10 = handles.spin - 1*handles.etxt_ff;          %
nspin20 = handles.spin - 2*handles.etxt_ff;          %
nspin30 = handles.spin - 3*handles.etxt_ff;          %

handles.ntheta10=ntheta10;
handles.ntheta20=ntheta20;
handles.ntheta30=ntheta30;
handles.nphi10=nphi10;
handles.nphi20=nphi20;
handles.nphi30=nphi30;
handles.nspin10=nspin10;
handles.nspin20=nspin20;
handles.nspin30=nspin30;        % lazy man's find/replace.
guidata(hObject,handles);       % update handles.

%--- negative zs by angle
%--- zs (theta)
nzx10t = sind(ntheta10)*cosd(handles.phi);           %
nzy10t = sind(ntheta10)*sind(handles.phi);           %
nzz10t = cosd(ntheta10);                             %
nzx20t = sind(ntheta20)*cosd(handles.phi);           %
nzy20t = sind(ntheta20)*sind(handles.phi);           %
nzz20t = cosd(ntheta20);                             %
nzx30t = sind(ntheta30)*cosd(handles.phi);           %
nzy30t = sind(ntheta30)*sind(handles.phi);           %
nzz30t = cosd(ntheta30);                             %
%--- zs (phi)
nzx10p = sind(handles.theta)*cosd(nphi10);           %
nzy10p = sind(handles.theta)*sind(nphi10);           %
nzz10p = cosd(handles.theta);                        %
nzx20p = sind(handles.theta)*cosd(nphi20);           %
nzy20p = sind(handles.theta)*sind(nphi20);           %
nzz20p = cosd(handles.theta);                        %
nzx30p = sind(handles.theta)*cosd(nphi30);           %
nzy30p = sind(handles.theta)*sind(nphi30);           %
nzz30p = cosd(handles.theta);                        %

%--- xs components by angle (theta)
nxx10t = sind(handles.ntheta10+90)*cosd(handles.phi+0);
nxy10t = sind(handles.ntheta10+90)*sind(handles.phi+0);
nxz10t = cosd(handles.ntheta10+90);
nxx20t = sind(handles.ntheta20+90)*cosd(handles.phi+0);
nxy20t = sind(handles.ntheta20+90)*sind(handles.phi+0);
nxz20t = cosd(handles.ntheta20+90);
nxx30t = sind(handles.ntheta30+90)*cosd(handles.phi+0);
nxy30t = sind(handles.ntheta30+90)*sind(handles.phi+0);
nxz30t = cosd(handles.ntheta30+90);
%--- xs component by angle (phi)
nxx10p = sind(handles.theta+90)*cosd(handles.nphi10+0);
nxy10p = sind(handles.theta+90)*sind(handles.nphi10+0);
nxz10p = cosd(handles.theta+90);
nxx20p = sind(handles.theta+90)*cosd(handles.nphi20+0);
nxy20p = sind(handles.theta+90)*sind(handles.nphi20+0);
nxz20p = cosd(handles.theta+90);
nxx30p = sind(handles.theta+90)*cosd(handles.nphi30+0);
nxy30p = sind(handles.theta+90)*sind(handles.nphi30+0);
nxz30p = cosd(handles.theta+90);
%--- xs (in-plane XY) component by angle (spin)
nxx10s = handles.xs*cosd(handles.nspin10)-handles.ys*sind(handles.nspin10);
nxy10s = handles.xs*sind(handles.nspin10)-handles.ys*sind(handles.nspin10);
% nxz10s = handles.xs(3);
nxx20s = handles.xs*cosd(handles.nspin20)-handles.ys*sind(handles.nspin20);
nxy20s = handles.xs*sind(handles.nspin20)-handles.ys*sind(handles.nspin20);
% nxz20s = handles.xs(3);
nxx30s = handles.xs*cosd(handles.nspin30)-handles.ys*sind(handles.nspin30);
nxy30s = handles.xs*sind(handles.nspin30)-handles.ys*sind(handles.nspin30);
% nxz30s = handles.xs(3);
guidata(hObject,handles);
disp(' ')                                            %
disp('--------- negative values --------')           %
disp(['theta used for calculating zs = ',num2str(handles.theta)]);
disp(['  phi used for calculating zs = ',num2str(handles.phi)]);
disp(['theta used for calculating gbg_nzs10 = ',num2str(ntheta10)]);
disp(['  phi used for calculating gbg_nzs10 = ',num2str(nphi10)]);
disp(['theta used for calculating gbg_nzs20 = ',num2str(ntheta20)]);
disp(['  phi used for calculating gbg_nzs20 = ',num2str(nphi20)]);
disp(['theta used for calculating gbg_nzs30 = ',num2str(ntheta30)]);
disp(['  phi used for calculating gbg_nzs30 = ',num2str(nphi30)]);
disp('--------end negative values------')            %
disp(' ')                                            %
% -10 theta zs (garbage)                             %
gbg_nzs10t = [nzx10t nzy10t nzz10t];                 %
gbg_nxs10t = [nxx10t nxy10t nxz10t];                 %
gbg_nys10t = cross(gbg_nzs10t,gbg_nxs10t);           %
% -20 theta zs (garbage)                             %
gbg_nzs20t = [nzx20t nzy20t nzz20t];                 %
gbg_nxs20t = [nxx20t nxy20t nxz20t];                 %
gbg_nys20t = cross(gbg_nzs20t,gbg_nxs20t);           %
% -30 theta zs (garbage)                             %
gbg_nzs30t = [nzx30t nzy30t nzz30t];                 %
gbg_nxs30t = [nxx10t nxy30t nxz30t];                 %
gbg_nys30t = cross(gbg_nzs30t,gbg_nxs30t);           %
% -10 phi zs (garbage)                               %
gbg_nzs10p = [nzx10p nzy10p nzz10p];                 %
gbg_nxs10p = [nxx10p nxy10p nxz10p];                 %
gbg_nys10p = cross(gbg_nzs10p,gbg_nxs10p);           %
% -20 phi zs (garbage)                               %
gbg_nzs20p = [nzx20p nzy20p nzz20p];                 %
gbg_nxs20p = [nxx20p nxy20p nxz20p];                 %
gbg_nys20p = cross(gbg_nzs20p,gbg_nxs20p);           %
% -30 phi zs (garbage)                               %
gbg_nzs30p = [nzx30p nzy30p nzz30p];                 %
gbg_nxs30p = [nxx30p nxy30p nxz30p];                 %
gbg_nys30p = cross(gbg_nzs30p,gbg_nxs30p);           %
%
% x & y vectors from spin                            %
handles.nbxs1=(handles.bc.*handles.nb)*nxx10s';      %
handles.nbxs2=(handles.bc.*handles.nb)*nxx20s';      %
handles.nbxs3=(handles.bc.*handles.nb)*nxx30s';      %
handles.nbys1=(handles.bc.*handles.nb)*nxy10s';      %
handles.nbys2=(handles.bc.*handles.nb)*nxy20s';      %
handles.nbys3=(handles.bc.*handles.nb)*nxy30s';      %
% negative x (nx) vectors                            %
handles.nbxt1 = handles.bc.*handles.nb*gbg_nxs10t';  %
handles.nbxt2 = handles.bc.*handles.nb*gbg_nxs20t';  %
handles.nbxt3 = handles.bc.*handles.nb*gbg_nxs30t';  %
handles.nbxp1 = handles.bc.*handles.nb*gbg_nxs10p';  %
handles.nbxp2 = handles.bc.*handles.nb*gbg_nxs20p';  %
handles.nbxp3 = handles.bc.*handles.nb*gbg_nxs30p';  %
% negative y (ny) vectors                            %
handles.nbyt1 = handles.bc.*handles.nb*gbg_nys10t';  %
handles.nbyt2 = handles.bc.*handles.nb*gbg_nys20t';  %
handles.nbyt3 = handles.bc.*handles.nb*gbg_nys30t';  %
handles.nbyp1 = handles.bc.*handles.nb*gbg_nys10p';  %
handles.nbyp2 = handles.bc.*handles.nb*gbg_nys20p';  %
handles.nbyp3 = handles.bc.*handles.nb*gbg_nys30p';  %
%-------END NEGATIVE VECTORS-------------------------%
guidata(hObject,handles);

