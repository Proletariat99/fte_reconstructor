function varargout = dveGUI_load_export(varargin)
% DVE_LOAD_DATA MATLAB code for dveGUI_load_export.fig
%      DVE_LOAD_DATA, by itself, creates a new DVE_LOAD_DATA or raises the existing
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
%      applied to the GUI before dveGUI_load_export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dveGUI_load_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dveGUI_load_export

% Last Modified by GUIDE v2.5 03-Aug-2011 12:14:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @dveGUI_load_export_OpeningFcn, ...
  'gui_OutputFcn',  @dveGUI_load_export_OutputFcn, ...
  'gui_LayoutFcn',  @dveGUI_load_export_LayoutFcn, ...
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


% --- Executes just before dveGUI_load_export is made visible.
function dveGUI_load_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dveGUI_load_export (see VARARGIN)

% Choose default command line output for dveGUI_load_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dveGUI_load_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = dveGUI_load_export_OutputFcn(hObject, eventdata, handles)
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


% --- Creates and returns a handle to the GUI figure. 
function h1 = dveGUI_load_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', 227.0068359375, ...
    'taginfo', struct(...
    'figure', 2, ...
    'uipanel', 2, ...
    'popupmenu', 4, ...
    'edit', 5, ...
    'text', 8, ...
    'pushbutton', 2), ...
    'override', 0, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', '/Users/dave/Documents/MATLAB/dveGUI_load_export.m', ...
    'lastFilename', '/Users/dave/Documents/MATLAB/dveGUI_load_data.fig');
appdata.lastValidTag = 'figure1';
appdata.GUIDELayoutEditor = [];
appdata.initTags = struct(...
    'handle', [], ...
    'tag', 'figure1');

h1 = figure(...
'Units','characters',...
'Color',[0.929411764705882 0.929411764705882 0.929411764705882],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','dveGUI_load_data',...
'NumberTitle','off',...
'PaperPositionMode','auto',...
'Position',[103.833333333333 36.8333333333333 25.6666666666667 24.6666666666667],...
'Resize','off',...
'HandleVisibility','callback',...
'UserData',[],...
'Tag','figure1',...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'data';

h2 = uipanel(...
'Parent',h1,...
'Units','characters',...
'FontUnits','pixels',...
'FontSize',10.0239122617763,...
'Title','Initial Data',...
'Tag','data',...
'UserData',[],...
'Clipping','on',...
'Position',[1.5 2.5 23 21.3333333333333],...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'pm_doy';

h3 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('pm_doy_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[6.33333333333333 17.6666666666667 15 2],...
'String',{  '2007/06/08'; '2007/11/14'; '2008/03/22'; '2008/04/05'; '2008/04/13'; '2008/04/15'; '2008/04/17'; '2008/04/19'; '2008/05/17'; '2008/11/30'; '2009/04/06'; '2009/04/23'; '2009/05/15'; '2009/05/29'; '2009/11/02' },...
'Style','popupmenu',...
'TooltipString','Format is 2010/05/08',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('pm_doy_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','pm_doy');

appdata = [];
appdata.lastValidTag = 'etxt_ti';

h4 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('etxt_ti_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[6.33333333333333 15 15 2],...
'String','06:30:00',...
'Style','edit',...
'TooltipString','Format is 16:34:56',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('etxt_ti_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','etxt_ti');

appdata = [];
appdata.lastValidTag = 'etxt_tf';

h5 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('etxt_tf_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[6.33333333333333 12.3333333333333 15 2],...
'String','08:30:00',...
'Style','edit',...
'TooltipString','Format is 17:34:56',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('etxt_tf_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','etxt_tf');

appdata = [];
appdata.lastValidTag = 'pm_sc';

h6 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('pm_sc_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[11.3333333333333 9.66666666666667 10 2],...
'String',{  'a'; 'b'; 'c' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('pm_sc_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','pm_sc');

appdata = [];
appdata.lastValidTag = 'pm_src';

h7 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('pm_src_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[5.83333333333333 7 15.5 2],...
'String',{  'reduced'; 'onboard' },...
'Style','popupmenu',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('pm_src_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','pm_src');

appdata = [];
appdata.lastValidTag = 'etxt_w';

h8 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('etxt_w_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[9.66666666666667 4.33333333333333 11.6666666666667 2],...
'String','2',...
'Style','edit',...
'TooltipString','Format is 10',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('etxt_w_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','etxt_w');

appdata = [];
appdata.lastValidTag = 'text1';

h9 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 18.0833333333333 5 1.08333333333333],...
'String','doy:',...
'Style','text',...
'UserData',[],...
'Tag','text1',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text2';

h10 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 12.75 5 1.08333333333333],...
'String','tf:',...
'Style','text',...
'UserData',[],...
'Tag','text2',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text3';

h11 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 7.41666666666667 5.33333333333333 1.5],...
'String','src:',...
'Style','text',...
'UserData',[],...
'Tag','text3',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text4';

h12 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 4.75 5 1.08333333333333],...
'String','w:',...
'Style','text',...
'UserData',[],...
'Tag','text4',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text5';

h13 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 10.0833333333333 5 1.08333333333333],...
'String','sc:',...
'Style','text',...
'UserData',[],...
'Tag','text5',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'text6';

h14 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 15.4166666666667 5 1.08333333333333],...
'String','ti:',...
'Style','text',...
'UserData',[],...
'Tag','text6',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'etxt_nev';

h15 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('etxt_nev_Callback',hObject,eventdata,guidata(hObject)),...
'CData',[],...
'FontSize',10.0239122617763,...
'Position',[14.5 1.66666666666667 6.83333333333333 2],...
'String','1',...
'Style','edit',...
'TooltipString','Format is 10',...
'CreateFcn', {@local_CreateFcn, @(hObject,eventdata)dve_load_data('etxt_nev_CreateFcn',hObject,eventdata,guidata(hObject)), appdata} ,...
'UserData',[],...
'Tag','etxt_nev');

appdata = [];
appdata.lastValidTag = 'text7';

h16 = uicontrol(...
'Parent',h2,...
'Units','characters',...
'FontUnits','pixels',...
'CData',[],...
'FontSize',11,...
'FontWeight','bold',...
'HorizontalAlignment','left',...
'Position',[0.5 2.08333333333333 12.6666666666667 1.08333333333333],...
'String','Event # (nev):',...
'Style','text',...
'UserData',[],...
'Tag','text7',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );

appdata = [];
appdata.lastValidTag = 'btn_load';

h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontUnits','pixels',...
'Callback',@(hObject,eventdata)dve_load_data('btn_load_Callback',hObject,eventdata,guidata(hObject)),...
'FontSize',10.0239122617763,...
'Position',[7.16666666666667 0.416666666666667 11.6666666666667 1.75],...
'String','Load Data',...
'Tag','btn_load',...
'CreateFcn', {@local_CreateFcn, blanks(0), appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   if isa(createfcn,'function_handle')
       createfcn(hObject, eventdata);
   else
       eval(createfcn);
   end
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)

gui_StateFields =  {'gui_Name'
    'gui_Singleton'
    'gui_OpeningFcn'
    'gui_OutputFcn'
    'gui_LayoutFcn'
    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('MATLAB:gui_mainfcn:FieldNotFound', 'Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % DVEGUI_LOAD_EXPORT
    % create the GUI only if we are not in the process of loading it
    % already
    gui_Create = true;
elseif local_isInvokeActiveXCallback(gui_State, varargin{:})
    % DVEGUI_LOAD_EXPORT(ACTIVEX,...)
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif local_isInvokeHGCallback(gui_State, varargin{:})
    % DVEGUI_LOAD_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = false;
else
    % DVEGUI_LOAD_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = true;
end

if ~gui_Create
    % In design time, we need to mark all components possibly created in
    % the coming callback evaluation as non-serializable. This way, they
    % will not be brought into GUIDE and not be saved in the figure file
    % when running/saving the GUI from GUIDE.
    designEval = false;
    if (numargin>1 && ishghandle(varargin{2}))
        fig = varargin{2};
        while ~isempty(fig) && ~ishghandle(fig,'figure')
            fig = get(fig,'parent');
        end
        
        designEval = isappdata(0,'CreatingGUIDEFigure') || isprop(fig,'__GUIDEFigure');
    end
        
    if designEval
        beforeChildren = findall(fig);
    end
    
    % evaluate the callback now
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else       
        feval(varargin{:});
    end
    
    % Set serializable of objects created in the above callback to off in
    % design time. Need to check whether figure handle is still valid in
    % case the figure is deleted during the callback dispatching.
    if designEval && ishghandle(fig)
        set(setdiff(findall(fig),beforeChildren), 'Serializable','off');
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end

    % Check user passing 'visible' P/V pair first so that its value can be
    % used by oepnfig to prevent flickering
    gui_Visible = 'auto';
    gui_VisibleInput = '';
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        % Recognize 'visible' P/V pair
        len1 = min(length('visible'),length(varargin{index}));
        len2 = min(length('off'),length(varargin{index+1}));
        if ischar(varargin{index+1}) && strncmpi(varargin{index},'visible',len1) && len2 > 1
            if strncmpi(varargin{index+1},'off',len2)
                gui_Visible = 'invisible';
                gui_VisibleInput = 'off';
            elseif strncmpi(varargin{index+1},'on',len2)
                gui_Visible = 'visible';
                gui_VisibleInput = 'on';
            end
        end
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.

    
    % Do feval on layout code in m-file if it exists
    gui_Exported = ~isempty(gui_State.gui_LayoutFcn);
    % this application data is used to indicate the running mode of a GUIDE
    % GUI to distinguish it from the design mode of the GUI in GUIDE. it is
    % only used by actxproxy at this time.   
    setappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]),1);
    if gui_Exported
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);

        % make figure invisible here so that the visibility of figure is
        % consistent in OpeningFcn in the exported GUI case
        if isempty(gui_VisibleInput)
            gui_VisibleInput = get(gui_hFigure,'Visible');
        end
        set(gui_hFigure,'Visible','off')

        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen');
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt, gui_Visible);
        end
    end
    if isappdata(0, genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]))
        rmappdata(0,genvarname(['OpenGuiWhenRunning_', gui_State.gui_Name]));
    end

    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    % Singleton setting in the GUI M-file takes priority if different
    gui_Options.singleton = gui_State.gui_Singleton;

    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA. If there is
        % user set GUI data already, keep that also.
        data = guidata(gui_hFigure);
        handles = guihandles(gui_hFigure);
        if ~isempty(handles)
            if isempty(data)
                data = handles;
            else
                names = fieldnames(handles);
                for k=1:length(names)
                    data.(char(names(k)))=handles.(char(names(k)));
                end
            end
        end
        guidata(gui_hFigure, data);
    end

    % Apply input P/V pairs other than 'visible'
    for index=1:2:length(varargin)
        if length(varargin) == index || ~ischar(varargin{index})
            break;
        end

        len1 = min(length('visible'),length(varargin{index}));
        if ~strncmpi(varargin{index},'visible',len1)
            try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
        end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end

    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        % Handle the default callbacks of predefined toolbar tools in this
        % GUI, if any
        guidemfile('restoreToolbarToolPredefinedCallback',gui_hFigure); 
        
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);

        % Call openfig again to pick up the saved visibility or apply the
        % one passed in from the P/V pairs
        if ~gui_Exported
            gui_hFigure = local_openfig(gui_State.gui_Name, 'reuse',gui_Visible);
        elseif ~isempty(gui_VisibleInput)
            set(gui_hFigure,'Visible',gui_VisibleInput);
        end
        if strcmpi(get(gui_hFigure, 'Visible'), 'on')
            figure(gui_hFigure);
            
            if gui_Options.singleton
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        if isappdata(gui_hFigure,'InGUIInitialization')
            rmappdata(gui_hFigure,'InGUIInitialization');
        end

        % If handle visibility is set to 'callback', turn it on until
        % finished with OutputFcn
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end

    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end

    if isscalar(gui_hFigure) && ishghandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end

function gui_hFigure = local_openfig(name, singleton, visible)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
if nargin('openfig') == 2
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
else
    gui_hFigure = openfig(name, singleton, visible);
end

function result = local_isInvokeActiveXCallback(gui_State, varargin)

try
    result = ispc && iscom(varargin{1}) ...
             && isequal(varargin{1},gcbo);
catch
    result = false;
end

function result = local_isInvokeHGCallback(gui_State, varargin)

try
    fhandle = functions(gui_State.gui_Callback);
    result = ~isempty(findstr(gui_State.gui_Name,fhandle.file)) || ...
             (ischar(varargin{1}) ...
             && isequal(ishghandle(varargin{2}), 1) ...
             && (~isempty(strfind(varargin{1},[get(varargin{2}, 'Tag'), '_'])) || ...
                ~isempty(strfind(varargin{1}, '_CreateFcn'))) );
catch
    result = false;
end


