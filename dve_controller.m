function varargout = dve_controller(varargin)
%DVE_CONTROLLER M-file for dve_controller.fig
%      DVE_CONTROLLER, by itself, creates a new DVE_CONTROLLER or raises the existing
%      singleton*.
%
%      H = DVE_CONTROLLER returns the handle to a new DVE_CONTROLLER or the handle to
%      the existing singleton*.
%
%      DVE_CONTROLLER('Property','Value',...) creates a new DVE_CONTROLLER using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to dve_controller_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DVE_CONTROLLER('CALLBACK') and DVE_CONTROLLER('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DVE_CONTROLLER.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_controller

% Last Modified by GUIDE v2.5 03-Aug-2011 15:04:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dve_controller_OpeningFcn, ...
                   'gui_OutputFcn',  @dve_controller_OutputFcn, ...
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


% --- Executes just before dve_controller is made visible.
function dve_controller_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for dve_controller
% handles.input = varargin;
handles.output = hObject;


% handles.input = varargin;
% Update handles structure
guidata(hObject, handles);
% UIWAIT makes dve_controller wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% disp(['handles.mvab_theta= ',num2str(handles.mvab_theta)])
% disp(['handles.mvab_phi= ',num2str(handles.mvab_phi)])



function varargout = dve_controller_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;



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
% handles.phi=str2double(get(hObject,'String'));
% disp(['phi from CreateFcn= ', num2str(handles.phi)])
% ldh = load('vars1.mat');
% handles.ldh=ldh;
% guidata(hObject,ldh)
% disp(['mvab_theta= ',num2str(handles.ldh.mvab_theta)]);
% handles.phi=str2double(handles.ldh.mvab_theta);
% set(hObject,'String',num2str(handles.phi));
% guidata(hObject, handles);


function etxt_theta_Callback(hObject, eventdata, handles)
% handles.theta=str2double(get(hObject,'String'));
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
  
% --- Executes on button press in btn_disp.
function btn_disp_Callback(hObject, eventdata, handles)
%   handles.phi=str2double(get(handles.etxt_phi,'String'));
%   handles.theta=str2double(get(handles.etxt_theta,'String'));
%   phi=guidata(handles.etxt_phi,'String')
%   handles.theta = guidata(handles.etxt_theta,'String')
%   nx=sind(handles.theta)*cosd(handles.phi);
%   ny=sind(handles.theta)*sind(handles.phi);
%   nz=cosd(handles.theta);
%   nitr=[nx,ny,nz];
%   handles.zs = nitr;
%   handles.ys = cross(zs,xs);
%   
%   save('vars2.mat')
%   dve_display(handles);
%   
%   guidata(hObject, handles);


function btn_break_Callback(hObject, eventdata, handles)
disp('BREAK')
