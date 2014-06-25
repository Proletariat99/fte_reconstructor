function varargout = gui_gs(varargin)
% GUI_GS MATLAB code for gui_gs.fig
%      GUI_GS, by itself, creates a new GUI_GS or raises the existing
%      singleton*.
%
%      H = GUI_GS returns the handle to a new GUI_GS or the handle to
%      the existing singleton*.
%
%      GUI_GS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_GS.M with the given input arguments.
%
%      GUI_GS('Property','Value',...) creates a new GUI_GS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_gs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_gs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_gs

% Last Modified by GUIDE v2.5 18-Jul-2011 10:53:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_gs_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_gs_OutputFcn, ...
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


% --- Executes just before gui_gs is made visible.
function gui_gs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_gs (see VARARGIN)

% Choose default command line output for gui_gs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_gs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_gs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_accept.
function pushbutton_accept_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_accept (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_reject.
function pushbutton_reject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
