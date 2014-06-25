function varargout = dve_GUIfinal_v1(varargin)
%DVE_GUIFINAL_V1 M-file for dve_GUIfinal_v1.fig
%      DVE_GUIFINAL_V1, by itself, creates a new DVE_GUIFINAL_V1 or raises the existing
%      singleton*.
%
%      H = DVE_GUIFINAL_V1 returns the handle to a new DVE_GUIFINAL_V1 or the handle to
%      the existing singleton*.
%
%      DVE_GUIFINAL_V1('Property','Value',...) creates a new DVE_GUIFINAL_V1 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to dve_GUIfinal_v1_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DVE_GUIFINAL_V1('CALLBACK') and DVE_GUIFINAL_V1('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVE_GUIFINAL_V1.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_GUIfinal_v1

% Last Modified by GUIDE v2.5 10-Aug-2011 19:58:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dve_GUIfinal_v1_OpeningFcn, ...
  'gui_OutputFcn',  @dve_GUIfinal_v1_OutputFcn, ...
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
function dve_GUIfinal_v1_OpeningFcn(hObject, eventdata, handles, varargin)
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
% handles = catstruct(handles,guihandles);
% disp(['Class of handles= ',class(handles)]);
% disp(['Size of handles= ',num2str(size(handles))]);


% Choose default command line output for dve_GUIfinal_v1
handles.output = hObject;

guidata(hObject, handles);
function varargout = dve_GUIfinal_v1_OutputFcn(hObject, eventdata, handles)
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
handles.msh_vx = handles.msh{13}(handles.nev);
handles.msh_vy = handles.msh{14}(handles.nev);
handles.msh_vz = handles.msh{15}(handles.nev);
handles.msh_v = [handles.msh_vx,handles.msh_vy,handles.msh_vz];
handles.V0=handles.msh_v;
handles.nV0=handles.V0/norm(handles.V0);
handles.xsv0 = -handles.nV0;
handles.secs = handles.a.data(:,1);
handles.BX=handles.a.data(:,2);
handles.BY=handles.a.data(:,3);
handles.BZ=handles.a.data(:,4); 	% Assigns Bx, By, Bz (from file - GSM)
handles.VX=handles.a.data(:,6);
handles.VY=handles.a.data(:,7);
handles.VZ=handles.a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
handles.Nx=handles.a.data(:,10); 					% assigns plasma N density from file (cm^-3)
handles.Tx=handles.a.data(:,11); 					% assigns temp from file (K)
idata=1:length(handles.BX);                 % index from 1 to (# of handles.BX)
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
handles.bc=[handles.BX(handles.i1:handles.i2) handles.BY(handles.i1:handles.i2) handles.BZ(handles.i1:handles.i2)];				% re-assigns B for some reason, should be the same as above.
handles.bx = handles.bc(:,1);
handles.by = handles.bc(:,2);
handles.bz = handles.bc(:,3);
handles.vp=[handles.VX(handles.i1:handles.i2) handles.VY(handles.i1:handles.i2) handles.VZ(handles.i1:handles.i2)];
handles.Nx = handles.Nx;
handles.Tx = handles.Tx;
handles.lenb = length(handles.bc(:,1));
% Density for window ta to tb
[EHT,EC,HTR,HTslope,vht]=HTcoef(handles.bc,handles.vp,handles.lenb);
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(handles.bc,handles.lenb);
% handles.xtemp = -Kbint;
% handles.ytemp = Kbmax;
% handles.ztemp = Kbmin;
handles.xtemp=Kbmin;
handles.ytemp=Kbmax;
handles.ztemp=-Kbint;
nxorig=handles.ztemp(1);
nyorig=handles.ztemp(2);
nzorig=handles.ztemp(3);
handles.zs=[nxorig nyorig nzorig];

handles.xstemp=handles.xsv0;
handles.zs=cross(handles.xs,handles.ystemp);
handles.ys=cross(handles.zs,handles.xs);
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
handles.phi=handles.mvab_phi;
set(handles.txt_phi,'String',num2str(handles.mvab_phi));
set(handles.etxt_phi,'String',num2str(handles.mvab_phi));
handles.theta=handles.mvab_theta;
set(handles.txt_theta,'String',num2str(handles.mvab_theta));
set(handles.etxt_theta,'String',num2str(handles.mvab_theta));
disp(['Theta=',num2str(handles.theta)])
disp(['Phi=',num2str(handles.phi)])
disp(['actual msh V0 (V0)=',num2str(handles.V0)])
disp(['normalized V0 (nV0)=',num2str(handles.nV0)])
disp(['xs from V0 (-nV0)=',num2str(handles.xs)])
disp(['ys from MVAB=',num2str(handles.ys)])
disp(['zs from cross=',num2str(handles.zs)])
disp(' ')
disp('--------- Orthogonal Test (should be very close to 0.0) ---------')
disp(['xs dot ys =', num2str(handles.xs*handles.ys')])
disp(['ys dot zs =', num2str(handles.ys*handles.zs')])
disp(['xs dot zs =', num2str(handles.xs*handles.zs')])
disp(' ')
disp('------ Spherical Coords of each Axis in GSM -----')
[xtheta,xphi,xr]=cart2sph(handles.xs(1),handles.xs(2),handles.xs(3));
[ytheta,yphi,yr]=cart2sph(handles.ys(1),handles.ys(2),handles.ys(3));
[ztheta,zphi,zr]=cart2sph(handles.zs(1),handles.zs(2),handles.zs(3));
disp(['xtheta= ',num2str(rad2deg(xtheta))]);
disp(['  xphi= ',num2str(rad2deg(xphi))]);
disp(['ytheta= ',num2str(rad2deg(ytheta))]);
disp(['  yphi= ',num2str(rad2deg(zphi))]);
disp(['ztheta= ',num2str(rad2deg(ztheta))]);
disp(['  zphi= ',num2str(rad2deg(zphi))]);
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
uicontrol(handles.etxt_azi);
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
uicontrol(handles.btm_display);
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
%--
% Control Panel
%--
function btn_prev_Callback(hObject, eventdata, handles)
function btn_next_Callback(hObject, eventdata, handles)
function btn_exit_Callback(hObject, eventdata, handles)
function btn_save_Callback(hObject, eventdata, handles)
fig = gcf;
        anm = num2str(azi);
        pnm = num2str(round(pol));
        dnm = num2str([yr,'_',month,'_',day]);
        thnm = num2str(round(theta));
        phnm = num2str(round(phi));
        tmnm = num2str([ta,'to',tb]);
        filename = [anm,'a',pnm,'p',thnm, 'th', phnm, 'ph', dnm,'_',tmnm, '.png'];
        saveas(fig,filename,'png');
        disp(['File saved as ', filename])
function btn_break_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
disp('remove break or ctrl-option-r to continue')

%----------- Display Button Callback -------------%
function btn_display_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla(handles.axes1);
guidata(hObject, handles);
handles.ndata = handles.lenb;
handles.dTim1 = handles.secs(handles.i1:handles.i2)-handles.secs(handles.i1);               % my code (dd, 2011)
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
% project B into reconstruction coord.
handles.bxs=(handles.bc.*handles.nb)*handles.xs';
handles.bys=(handles.bc.*handles.nb)*handles.ys';
handles.bzs=(handles.bc.*handles.nb)*handles.zs';
handles.vxs = -handles.V0*handles.xs';                  % re-adjust V0 into new xs
handles.xa0=handles.vxs.*handles.dTim1;                 % calculate distance (V*t)
handles.iix = handles.i1:handles.i2;                    % index of interesting values.
handles.dzn=handles.Nx(handles.iix);                    % density reconstruction by matrix
handles.Tp = handles.Tx(handles.iix)/handles.kb2;       % temperature selection and mod
%Project V into reconstruction coord.
vxc=zeros(1,handles.ndata);                             % zeros (for speed)
vyc=zeros(1,handles.ndata);
vzc=zeros(1,handles.ndata);
for i=1:handles.ndata
  vxc(i)=(handles.vp(i,:)-handles.V0)*handles.xs';      % projects V into new x, y & z
  vyc(i)=(handles.vp(i,:)-handles.V0)*handles.ys';
  vzc(i)=(handles.vp(i,:)-handles.V0)*handles.zs';
end
%Calculate vector potential A along y=0

%-------------------%
% Note: xs = old zs %
%       ys = old xs %
%       zs = old ys %
%-------------------%
ndA=zeros(1,handles.lenb);                              
for i=2:handles.lenb                       %
%   ndA(i)=-(handles.bys(i)+handles.bys(i-1))*(handles.dTim1(i)-handles.dTim1(i-1))*handles.vxs*handles.nv*0.5;  % Determines A for different
    ndA(i)=-(handles.bzs(i)+handles.bzs(i-1))*(handles.dTim1(i)-handles.dTim1(i-1))*handles.vxs*handles.nv*0.5;  % Determin
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
% handles.pbz=pp./p0+((handles.bzs./handles.b0).^2)/2;
handles.pbx=pp./p0+((handles.bxs./handles.b0).^2)/2;
handles.An=A1./A0;
guidata(hObject, handles);
handles.amax1=max(handles.An);
handles.amin1=min(handles.An);
handles.fS1=polyfit(handles.An,handles.pbx,4); %Pt(A)
% handles.fZ1=polyfit(handles.An,handles.bzs/handles.b0,5); %Bz(A)
handles.fX1=polyfit(handles.An,handles.bxs/handles.b0,5); %Bz(A)
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
handles.axes1;
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
cla(handles.axes6);
btn_vectorload_Callback(hObject,eventdata,handles)
phick=get(handles.cb_showphi,'Value');
thtck=get(handles.cb_showtheta,'Value');
ckarr=[phick thtck];
str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)));
switch str_ck
  case '00'
    cla(handles.axes6);
  case '11'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showphi_Callback(hObject,eventdata,handles);
    cb_showtheta_Callback(hObject,eventdata,handles);
  case '10'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showphi_Callback(hObject,eventdata,handles);
  case '01'
    cla(handles.axes6);
    btn_vectorload_Callback(hObject,eventdata,handles)
    cb_showtheta_Callback(hObject,eventdata,handles);
  otherwise
    disp('Some other issue exists. Oh, and go brush up your linear algebra, Captain Clunky.')
end
disp(' ')
disp('oooooooooooooooooo-- End Display Button Data --oooooooooooooooooo')
function plotstuff(hObject, eventdata,handles)
[Aup,Aupdy,Aupdx]=gsup(handles.c1,handles.x,handles.y,handles.u,handles.udy,handles.udx,handles.fS1);
[Adn,Adndy,Adndx]=gsdn(handles.c1,handles.x,handles.y,handles.u,handles.udy,handles.udx,handles.fS1);
Aup=Aup+Adn;
Aup(handles.mid,:)=Aup(handles.mid,:)/2;
Aupdy=Aupdy+Adndy;
Aupdy(handles.mid,:)=Aupdy(handles.mid,:)/2;
Aupdx=Aupdx+Adndx;
Aupdx(handles.mid,:)=Aupdx(handles.mid,:)/2;
handles.gg=1:handles.ny;
Bzup=zeros(handles.ny,handles.nx);
for j=1:handles.ny
  for i=1:handles.nx
    Bzup(j,i)=polyval(handles.fX1,Aup(j,i));
    if Aup(j,i)>handles.amax1
      Bzup(j,i)=polyval(handles.fX1,handles.amax1);
    elseif Aup(j,i)<handles.amin1
      Bzup(j,i)=polyval(handles.fX1,handles.amin1);
    end
  end
end
B3 = Bzup(handles.gg,:)*handles.b0*1e9;              % defines 3rd axis
pcolor(handles.x*handles.L0,handles.y(handles.gg)*handles.L0,B3);            % displays 3rd axis (bz) in faceted shading
shading interp                       % changes shading to interpolated
minb3=min(min(B3));
maxb3=max(max(B3));
caxis([minb3 maxb3]);
zbar=colorbar('vertical');
set(zbar,'Fontsize',8);
set(get(zbar,'XLabel'),'String','Bz [nT]','Rotation',0,'HorizontalAlignment','left','Fontsize',10);
hold on;
[cc,hh]=contour(handles.x*handles.L0,handles.y(handles.gg)*handles.L0,Aup(handles.gg,:),[-4:0.11:4],'k');
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
  btn_vectorload_Callback(hObject,eventdata,handles)
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
  cla(handles.axes6);
  btn_vectorload_Callback(hObject,eventdata,handles)
  phick=get(handles.cb_showphi,'Value');
  thtck=get(handles.cb_showtheta,'Value');
  ckarr=[phick thtck];
  str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)));
  switch str_ck
    case '00'
      cla(handles.axes6);
    case '11'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
    case '10'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
    case '01'
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
  btn_vectorload_Callback(hObject,eventdata,handles)
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
  btn_vectorload_Callback(hObject,eventdata,handles)
  phick=get(handles.cb_showphi,'Value');
  thtck=get(handles.cb_showtheta,'Value');
  ckarr=[phick thtck];
  str_ck=strcat(num2str(ckarr(1)), num2str(ckarr(2)));
  switch str_ck
    case '00'
      cla(handles.axes6);
    case '11'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
      cb_showtheta_Callback(hObject,eventdata,handles);
    case '10'
      cla(handles.axes6);
      btn_vectorload_Callback(hObject,eventdata,handles)
      cb_showphi_Callback(hObject,eventdata,handles);
    case '01'
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
%--- xs (formerly zs) components by angle
%--- +10 theta
xx10t = sind(handles.theta10)*cosd(handles.phi);
xy10t = sind(handles.theta10)*sind(handles.phi);
xz10t = cosd(handles.theta10);
%--- +20 theta
xx20t = sind(handles.theta20)*cosd(handles.phi);
xy20t = sind(handles.theta20)*sind(handles.phi);
xz20t = cosd(handles.theta20);
%--- +30 theta
xx30t = sind(handles.theta30)*cosd(handles.phi);
xy30t = sind(handles.theta30)*sind(handles.phi);
xz30t = cosd(handles.theta30);
%--- +10 phi
xx10p = sind(handles.theta)*cosd(handles.phi20);
xy10p = sind(handles.theta)*sind(handles.phi20);
xz10p = cosd(handles.theta);
%--- +20 phi
xx20p = sind(handles.theta)*cosd(handles.phi20);
xy20p = sind(handles.theta)*sind(handles.phi20);
xz20p = cosd(handles.theta);
%--- +30 phi
xx30p = sind(handles.theta)*cosd(handles.phi30);
xy30p = sind(handles.theta)*sind(handles.phi30);
xz30p = cosd(handles.theta);
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
%--- ys (formerly xs) components by angle (theta)
yx10t = cosd(handles.theta10)*cosd(handles.phi);
yy10t = cosd(handles.theta10)*sind(handles.phi);
yz10t = -sind(handles.theta10);
yx20t = cosd(handles.theta20)*cosd(handles.phi);
yy20t = cosd(handles.theta20)*sind(handles.phi);
yz20t = -sind(handles.theta20);
yx30t = cosd(handles.theta30)*cosd(handles.phi);
yy30t = cosd(handles.theta30)*sind(handles.phi);
yz30t = -sind(handles.theta30);
%--- ys (formerly xs) components by angle (phi)
yx10p = cosd(handles.theta)*cosd(handles.phi10);
yy10p = cosd(handles.theta)*sind(handles.phi10);
yz10p = -sind(handles.theta);
yx20p = cosd(handles.theta)*cosd(handles.phi20);
yy20p = cosd(handles.theta)*sind(handles.phi20);
yz20p = -sind(handles.theta);
yx30p = cosd(handles.theta)*cosd(handles.phi30);
yy30p = cosd(handles.theta)*sind(handles.phi30);
yz30p = -sind(handles.theta);
%--- zs by angle (do we need this?)
% +10 theta xs (gbg = garbage)
gbg_xs10t = [xx10t xy10t xz10t];
gbg_ys10t = [yx10t yy10t yz10t];
gbg_zs10t = cross(gbg_xs10t,gbg_ys10t);
% +20 theta zs (garbage)
gbg_xs20t = [xx20t xy20t xz20t];
gbg_ys20t = [yx20t yy20t yz20t];
gbg_zs20t = cross(gbg_xs20t,gbg_ys20t);
% + 30 theta zs (garbage)
gbg_xs30t = [xx30t xy30t xz30t];
gbg_ys30t = [yx30t yy30t yz30t];
gbg_zs30t = cross(gbg_xs30t,gbg_ys30t);
% +10 phi zs (garbage)
gbg_xs10p = [xx10p xy10p xz10p];
gbg_ys10p = [yx10p yy10p yz10p];
gbg_zs10p = cross(gbg_xs10p,gbg_ys10p);
% +20 phi zs (garbage)
gbg_xs20p = [xx20p xy20p xz20p];
gbg_ys20p = [yx20p yy20p yz20p];
gbg_zs20p = cross(gbg_xs20p,gbg_ys20p);
% + 30 phi zs (garbage)
gbg_xs30p = [xx30p xy30p xz30p];
gbg_ys30p = [yx30p yy30p yz30p];
gbg_zs30p = cross(gbg_xs30p,gbg_ys30p);
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

handles.ntheta10=ntheta10;
handles.ntheta20=ntheta20;
handles.ntheta30=ntheta30;
handles.nphi10=nphi10;
handles.nphi20=nphi20;
handles.nphi30=nphi30;      % lazy man's find/replace%
guidata(hObject,handles);       % update handles.

%--- negative xs by angle
%--- xs (formerly zs) (theta)
nxx10t = sind(ntheta10)*cosd(handles.phi);           %
nxy10t = sind(ntheta10)*sind(handles.phi);           %
nxz10t = cosd(ntheta10);                             %
nxx20t = sind(ntheta20)*cosd(handles.phi);           %
nxy20t = sind(ntheta20)*sind(handles.phi);           %
nxz20t = cosd(ntheta20);                             %
nxx30t = sind(ntheta30)*cosd(handles.phi);           %
nxy30t = sind(ntheta30)*sind(handles.phi);           %
nxz30t = cosd(ntheta30);                             %
%--- xs (formerly zs) (phi)
nxx10p = sind(handles.theta)*cosd(nphi10);           %
nxy10p = sind(handles.theta)*sind(nphi10);           %
nxz10p = cosd(handles.theta);                        %
nxx20p = sind(handles.theta)*cosd(nphi20);           %
nxy20p = sind(handles.theta)*sind(nphi20);           %
nxz20p = cosd(handles.theta);                        %
nxx30p = sind(handles.theta)*cosd(nphi30);           %
nxy30p = sind(handles.theta)*sind(nphi30);           %
nxz30p = cosd(handles.theta);                        %

%--- ys (formerly xs) components by angle (theta)
nyx10t = cosd(handles.ntheta10)*cosd(handles.phi);
nyy10t = cosd(handles.ntheta10)*sind(handles.phi);
nyz10t = -sind(handles.ntheta10);
nyx20t = cosd(handles.ntheta20)*cosd(handles.phi);
nyy20t = cosd(handles.ntheta20)*sind(handles.phi);
nyz20t = -sind(handles.ntheta20);
nyx30t = cosd(handles.ntheta30)*cosd(handles.phi);
nyy30t = cosd(handles.ntheta30)*sind(handles.phi);
nyz30t = -sind(handles.ntheta30);
%--- ys (formerly xs) component by angle (nphi)
nyx10p = cosd(handles.theta)*cosd(handles.nphi10);
nyy10p = cosd(handles.theta)*sind(handles.nphi10);
nyz10p = -sind(handles.theta);
nyx20p = cosd(handles.theta)*cosd(handles.nphi20);
nyy20p = cosd(handles.theta)*sind(handles.nphi20);
nyz20p = -sind(handles.theta);
nyx30p = cosd(handles.theta)*cosd(handles.nphi30);
nyy30p = cosd(handles.theta)*sind(handles.nphi30);
nyz30p = -sind(handles.theta);
guidata(hObject,handles);
disp(' ')                                            %
disp('--------- negative values --------')           %
disp(['theta used for calculating zs = ',num2str(handles.theta)]);
disp(['  phi used for calculating zs = ',num2str(handles.phi)]);
disp(['theta used for calculating gbg_nxs10 = ',num2str(ntheta10)]);
disp(['  phi used for calculating gbg_nxs10 = ',num2str(nphi10)]);
disp(['theta used for calculating gbg_nxs20 = ',num2str(ntheta20)]);
disp(['  phi used for calculating gbg_nxs20 = ',num2str(nphi20)]);
disp(['theta used for calculating gbg_nxs30 = ',num2str(ntheta30)]);
disp(['  phi used for calculating gbg_nxs30 = ',num2str(nphi30)]);
disp('--------end negative values------')            %
disp(' ')                                            %
% -10 theta zs (garbage)                             %
gbg_nys10t = [nyx10t nyy10t nyz10t];                 %
gbg_nxs10t = [nxx10t nxy10t nxz10t];                 %
gbg_nzs10t = cross(gbg_nxs10t,gbg_nys10t);           %
% -20 theta zs (garbage)                             %
gbg_nys20t = [nyx20t nyy20t nyz20t];                 %
gbg_nxs20t = [nxx20t nxy20t nxz20t];                 %
gbg_nzs20t = cross(gbg_nxs20t,gbg_nys20t);            %
% -30 theta zs (garbage)                             %
gbg_nys30t = [nyx30t nyy30t nyz30t];                 %
gbg_nxs30t = [nxx30t nxy30t nxz30t];                 %
gbg_nzs30t = cross(gbg_nxs30t,gbg_nys30t);            %
% -10 phi zs (garbage)                               %
gbg_nys10p = [nyx10p nyy10p nyz10p];                 %
gbg_nxs10p = [nxx10p nxy10p nxz10p];                 %
gbg_nzs10p = cross(gbg_nxs10p,gbg_nys10p);           %
% -20 phi zs (garbage)                               %
gbg_nys20p = [nyx20p nyy20p nyz20p];                 %
gbg_nxs20p = [nxx20p nxy20p nxz20p];                 %
gbg_nzs20p = cross(gbg_nxs20p,gbg_nys20p);          %
% -30 phi zs (garbage)                               %
gbg_nys30p = [nyx30p nyy30p nyz30p];                 %
gbg_nxs30p = [nxx30p nxy30p nxz30p];                 %
gbg_nzs30p = cross(gbg_nxs30p,gbg_nys30p);           %
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
