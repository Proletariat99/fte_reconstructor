function varargout = dve_load_data(varargin)
%DVE_LOAD_DATA M-file for dve_load_data.fig
%      DVE_LOAD_DATA, by itself, creates a new DVE_LOAD_DATA or raises the existing
%      singleton*.
%
%      H = DVE_LOAD_DATA returns the handle to a new DVE_LOAD_DATA or the handle to
%      the existing singleton*.
%
%      DVE_LOAD_DATA('Property','Value',...) creates a new DVE_LOAD_DATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to dve_load_data_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DVE_LOAD_DATA('CALLBACK') and DVE_LOAD_DATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVE_LOAD_DATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_load_data

% Last Modified by GUIDE v2.5 04-Aug-2011 11:25:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dve_load_data_OpeningFcn, ...
                   'gui_OutputFcn',  @dve_load_data_OutputFcn, ...
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


% --- Executes just before dve_load_data is made visible.
function dve_load_data_OpeningFcn(hObject, eventdata, handles, varargin)
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

% Choose default command line output for dve_load_data
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
%movegui(dve_load_data,'northwest')
% UIWAIT makes dve_load_data wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dve_load_data_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%====================================%
% INTERFACE CALLBACK AND CREATE CODE %
%====================================%


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
guidata(hObject, handles);

function etxt_ti_Callback(hObject, eventdata, handles)
handles.ti=get(hObject,'String');
disp(['handles.ti from callbackFcn= ',handles.ti])
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
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

% --- Executes during object creation, after setting all properties.
function etxt_tf_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
end
handles.tf=get(hObject,'String');
disp(['handles.tf from createFcn= ',handles.tf])
guidata(hObject, handles);

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

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<%
%      LOAD DATA BUTTON CALLBACK       %
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>%
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
    handles.bV0=msh_v;
    handles.V0=handles.bV0/norm(handles.bV0);
    handles.xsv0 = -handles.V0;
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
    handles.vp=[VX(handles.i1:handles.i2) VY(handles.i1:handles.i2) VZ(handles.i1:handles.i2)];
    handles.Nx = handles.Nx;
    handles.Tx = handles.Tx;
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
    handles.zs = Kbmin;
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
    guidata(hObject, handles);


function etxt_phi_Callback(hObject, eventdata, handles)
handles.phi=str2double(get(hObject,'String'));
disp(['phi from CallbackFcn= ', num2str(handles.phi)])
guidata(hObject, handles);
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
guidata(hObject, handles);          %update again


function etxt_phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.phi=str2double(get(hObject,'String'));
disp(['phi from CreateFcn= ', num2str(handles.phi)])
% guidata(hObject,ldh)
% disp(['mvab_theta= ',num2str(handles.ldh.mvab_theta)]);


function etxt_theta_Callback(hObject, eventdata, handles)
handles.theta=str2double(get(hObject,'String'));
disp(['theta from CallbackFcn= ', num2str(handles.theta)])
guidata(hObject, handles);          %update handles.
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
guidata(hObject, handles);          %update again


% --- Executes during object creation, after setting all properties.
function etxt_theta_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% handles.theta=str2double(get(hObject,'String'));
% disp(['handles.theta from CreateFcn= ', num2str(handles.theta)])
% disp(['handles.mvab_theta= ',num2str(handles.mvab_theta)])
guidata(hObject, handles);

function etxt_azi_Callback(hObject, eventdata, handles)
% hObject    handle to etxt_azi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of etxt_azi as text
%        str2double(get(hObject,'String')) returns contents of etxt_azi as a double
handles.azi=str2double(get(hObject,'String'));
disp(['handles.azi from CallbackFcn= ', num2str(handles.azi)])
guidata(hObject, handles);
set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
guidata(hObject, handles);          %update again

% --- Executes during object creation, after setting all properties.
function etxt_azi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etxt_azi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.azi=str2double(get(hObject,'String'));
disp(['azi from CreateFcn= ', num2str(handles.azi)])

guidata(hObject, handles);

function etxt_pol_Callback(hObject, eventdata, handles)
handles.pol=str2double(get(hObject,'String'));
disp(['handles.pol from CallbackFcn= ', num2str(handles.pol)])
guidata(hObject, handles);
set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
guidata(hObject, handles);          %update again

% --- Executes during object creation, after setting all properties.
function etxt_pol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to etxt_pol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
handles.pol=str2double(get(hObject,'String'));
disp(['handles.pol from CreateFcn= ', num2str(handles.pol)])

guidata(hObject, handles);

% --- Executes on button press in btn_prev.
function btn_prev_Callback(hObject, eventdata, handles)
% hObject    handle to btn_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_next.
function btn_next_Callback(hObject, eventdata, handles)
% hObject    handle to btn_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_exit.
function btn_exit_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btn_phiup.
function btn_phiup_Callback(hObject, eventdata, handles)
  handles.phi=handles.phi+handles.azi;
  set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
  guidata(hObject, handles);
  
% --- Executes on button press in btn_phidn.
function btn_phidn_Callback(hObject, eventdata, handles)
  handles.phi=handles.phi-handles.azi;
  set(uicontrol(handles.txt_phi),'String',num2str(handles.phi));
  guidata(hObject, handles);

% --- Executes on button press in btn_thetaup.
function btn_thetaup_Callback(hObject, eventdata, handles)
  handles.theta=handles.theta+handles.pol;
  set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
  guidata(hObject, handles);

% --- Executes on button press in btn_thetadn.
function btn_thetadn_Callback(hObject, eventdata, handles)
  handles.theta=handles.theta-handles.pol;
  set(uicontrol(handles.txt_theta),'String',num2str(handles.theta));
  guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function txt_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called  

%---------------------------------------------%
  %            GRAPH BUTTON CALLBACK            %
    %---------------------------------------------%
function btn_graph_Callback(hObject, eventdata, handles)
display_graph();

%===================================================================%
%                      DISPLAY BUTTON CALLBACK                      %
%===================================================================%
function btn_display_Callback(hObject, eventdata, handles)
axes(handles.axes1);
cla;
% Importing Handles variables to local variables (clunky, but fast);
guidata(hObject, handles);
bc = handles.bc;
nb = handles.nb;
nv = handles.nv;
kb = handles.kb;
kb2 = handles.kb2;
nden = handles.nden;
miu = handles.miu;
mp = handles.mp;
ndata = handles.lenb;
vp = handles.vp;
i1 = handles.i1;
i2 = handles.i2;
secs = handles.secs;
handles.dTim1 = secs(i1:i2)-secs(i1);               % my code (dd, 2011)
dTim1 = handles.dTim1;
V0 = handles.V0;
len = handles.lenb;
zx = sind(handles.theta)*cosd(handles.phi);
zy = sind(handles.theta)*sind(handles.phi);
zz = cosd(handles.theta);
handles.zs = [zx zy zz];
handles.ys = cross(handles.zs,handles.xsv0);
handles.xs = cross(handles.ys, handles.zs);     % recalculate x so all axes are perp.
xs = handles.xs;
ys = handles.ys;
zs = handles.zs;
Nx = handles.Nx;
bxs=(bc.*nb)*xs';
bys=(bc.*nb)*ys';
bzs=(bc.*nb)*zs';
vxs = -V0*xs';
xa0=vxs.*dTim1;
handles.iix = i1:i2;
iix = handles.iix;
handles.dzn=Nx(iix);                    % density reconstruction by matrix
dzn = handles.dzn;
Tx = handles.Tx;
Nx = handles.Nx;
handles.Tp = Tx(iix)/kb2;
Tp = handles.Tp;


%Project V into reconstruction coord.
vxc=zeros(1,ndata);
vyc=zeros(1,ndata);
vzc=zeros(1,ndata);
for i=1:ndata
  vxc(i)=(vp(i,:)-V0)*xs';
  vyc(i)=(vp(i,:)-V0)*ys';
  vzc(i)=(vp(i,:)-V0)*zs';
end

%Calculate vector potential A along y=0

ndA=zeros(1,len);
for i=2:len                       %
  ndA(i)=-(bys(i)+bys(i-1))*(dTim1(i)-dTim1(i-1))*vxs*nv*0.5;  % Determines A for different 
end
A1=cumsum(ndA);                   % with every additional deltat, we add more area
A1=A1';                           % transpose A1
bb=sqrt(bxs.^2+bys.^2+bzs.^2);    % magnitude of in-plane bfield vector
pp=dzn.*Tp*kb*nden;                     % new pp factor (without nt)
bmax=max(bb);                           % maximum b magnitude
nmax=max(dzn)*nden;                     % density maximum value
%Normalized factors
b0=bmax;                                % max b arrow
p0=bmax*bmax/miu;                       % max pressure value
n0=nmax;                                % max density value from 5 lines above.
% T0=p0/(kb*n0);                        % unused
vv0=sqrt(b0^2/(miu*n0*mp)*1e-6);
A0=max(abs(A1));
L0=(A0/bmax)*1e-3;                      % L0 is used for everything.  What is it?
pbz=pp./p0+((bzs./b0).^2)/2;
An=A1./A0;
  guidata(hObject, handles);
%   display_bfield(handles);
% 
%     %-------------------------------%
%   %     display bfield function   %       Can't get to work as a sep function
% %-------------------------------%         Just using inline with disp callback 
% function display_bfield(handles)          % instead.
%guidata(hObject, handles);

amax1=max(An); 
amin1=min(An);
fS1=polyfit(An,pbz,4); %Pt(A)
fZ1=polyfit(An,bzs/b0,5); %Bz(A)
  %--interpolation--%
  nx=ndata+(ndata-1)*3;
  xi=xa0(1):(xa0(ndata)-xa0(1))/(nx-1):xa0(ndata);
  
  bxi=interp1(xa0,bxs./b0,xi,'spline');
  byi=interp1(xa0,bys./b0,xi,'spline');
  % bzi=interp1(xa0,bzs./b0,xi,'spline'); % unused here.
  
  %--------------------------%
  ny=151;
  %--------------------------%  wtf is ny and why is it hard-coded?
  
  ht=0;
  py=0.1/1;
  mid=round(ny/2)+ht;
  
  x=xi./L0;
  hx=x(2)-x(1);%uniform grids
  hy=py*hx;
  y = zeros(1,ny);
  for j=1:ny
    y(j)=(j-mid)*hy;
  end
  
  %---------------------------%
  clear('ndA');
  ndA(1)=0;
  for i=2:nx
    ndA(i)=-(byi(i)+byi(i-1))*(x(i)-x(i-1))*0.5;
  end
  A2=cumsum(ndA);
  u=zeros(ny,nx);
  udy=zeros(ny,nx);
  udx=zeros(ny,nx);
  u(mid,:)=A2;
  udy(mid,:)=bxi;
  udx(mid,:)=byi;  
  c1=[nx ny mid hx hy amax1 amin1];
  
  %----------------------------------------------------------------------------%
  %=============================== SubBlock 1 =================================%
  %----------------------------------------------------------------------------%
   
    [Aup,Aupdy,Aupdx]=gsup(c1,x,y,u,udy,udx,fS1);
    [Adn,Adndy,Adndx]=gsdn(c1,x,y,u,udy,udx,fS1);
    
    Aup=Aup+Adn;
    Aup(mid,:)=Aup(mid,:)/2;
    
    Aupdy=Aupdy+Adndy;
    Aupdy(mid,:)=Aupdy(mid,:)/2;
    
    Aupdx=Aupdx+Adndx;
    Aupdx(mid,:)=Aupdx(mid,:)/2;
    
    gg=1:ny;                
    Bzup=zeros(ny,nx);      
    for j=1:ny              
      for i=1:nx
        Bzup(j,i)=polyval(fZ1,Aup(j,i));
        if Aup(j,i)>amax1
          Bzup(j,i)=polyval(fZ1,amax1);
        elseif Aup(j,i)<amin1
          Bzup(j,i)=polyval(fZ1,amin1);
        end
      end
    end
    
    B3 = Bzup(gg,:)*b0*1e9;              % defines 3rd axis
    pcolor(x*L0,y(gg)*L0,B3);            % displays 3rd axis (bz) in faceted shading
    shading interp                       % changes shading to interpolated

    
    minb3=min(min(B3));
    maxb3=max(max(B3));
    caxis([minb3 maxb3]);
    
    hbar=colorbar('vertical');
    pos=get(hbar,'Position');
    set(hbar,'Position',[pos(1)+0.1 pos(2)+0.1 0.02 0.2],'Fontsize',8)
    set(get(hbar,'XLabel'),'String','Bz [nT]','Rotation',0,...
      'HorizontalAlignment','left','Fontsize',8);
    
     hold on
    [cc,hh]=contour(x*L0,y(gg)*L0,Aup(gg,:),[-4:0.11:4],'k');
    set(hh,'linewidth',1.0);
    
    quiv1=quiver(xa0,zeros(1,ndata),bxs',bys',0.3,'w');  % plots white arrows accross middle.
    set(quiv1,'linewidth',1.0);                          % sets linewidth for quiver arrows
    
    set(handles.axes1,'fontsize',10,'TickDir','out','linewidth',1.0);
    axis equal
    xdist = x*L0;
    if xdist <= 0
        axis([min(x*L0), 0, min(y*L0), max(y*L0)])
    else
        axis([0, max(x*L0), min(y*L0), max(y*L0)])
     end
    
    % axis([0 500 min(y*L0) max(y*L0)])
    xlabel('x [km]','fontsize',14)
    ylabel('y [km]','fontsize',14)


    %-------------------------------%
  %     display graphs function   %
%-------------------------------%
% function display_graphs(hObject,eventdata,handles)
  axes(handles.axes2)
  cla;
  xm=1:ndata;
  plot(xm,vxc./vv0,'k-o',xm,vyc./vv0,'b-o',xm,vzc./vv0,'r-o','Markersize',3); grid on;
  ylabel('V','rotation',0);
  set(handles.axes2,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
  %--
 
  axes(handles.axes4)
  cla;
  plot(xm,bxs./b0,'k-o',xm,bys./b0,'b-o',xm,bzs./b0,'r-o','Markersize',3); grid on;
  ylabel('B','rotation',0);
  set(handles.axes4,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
  %--
  axes(handles.axes5)
  cla;
  plot(An,'k-o','Markersize',3);  grid on;
  ylabel('A','rotation',0);
  set(handles.axes5,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
  %--
  %Curve Fitting Processing
  %--
  fS1=polyfit(An,pbz,4); %Pt(A)
  fZ1=polyfit(An,bzs/b0,5); %Bz(A)
  
  amax1=max(An); amin1=min(An);
  
  ax1=amin1:0.01:amax1;
  
  df1=polyder(fS1); %dPt/dA
  yys1=polyval(fS1,ax1);
  
  yyz1=polyval(fZ1,ax1);

  axes(handles.axes6)
  cla;
  plot(An,pbz,'ko','Markersize',5); grid on;
   hold on;
  plot(ax1,yys1,'color',[0.5 0.5 0.5],'linewidth',3);
  set(handles.axes6,'ylim',[min(pbz)-0.1 max(pbz)+0.1],'ytick',[0:0.1:1],'fontsize',8);
  ylabel('p_t','rotation',0);
  xlabel('A')
  
  axes(handles.axes7)
  cla;
  plot(ax1,polyval(df1,ax1),'k'); grid on;
  ylabel('Dpbz');
  set(handles.axes7,'fontsize',8);
  
  axes(handles.axes8)
  cla;
  plot(An,bzs/b0,'ko','Markersize',3); grid on;
   hold on;
  plot(ax1,yyz1,'b','linewidth',3);
  set(handles.axes8,'ylim',[min(bzs/b0)-0.1 max(bzs/b0)+0.1],'ytick',[-1:0.25:1],'fontsize',8);
  ylabel('Bz');
  xlabel('A')


function btn_break_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
disp('keep this break in place')
