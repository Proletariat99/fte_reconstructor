function varargout = dve_display(varargin)
% DVE_DISPLAY MATLAB code for dve_display.fig
%      DVE_DISPLAY, by itself, creates a new DVE_DISPLAY or raises the existing
%      singleton*.
%
%      H = DVE_DISPLAY returns the handle to a new DVE_DISPLAY or the handle to
%      the existing singleton*.
%
%      DVE_DISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DVE_DISPLAY.M with the given input arguments.
%
%      DVE_DISPLAY('Property','Value',...) creates a new DVE_DISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dve_display_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dve_display_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dve_display

% Last Modified by GUIDE v2.5 03-Aug-2011 12:38:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dve_display_OpeningFcn, ...
                   'gui_OutputFcn',  @dve_display_OutputFcn, ...
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

% --- Executes just before dve_display is made visible.
function dve_display_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for dve_display
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using dve_display.
if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

% UIWAIT makes dve_display wait for user response (see UIRESUME)
% uiwait(handles.figure1);
ndata = handles.lenb;



%=======================================%
% Begin Original GS Reconstruction Code %
%=======================================%
%--interpolation--%
   nx=ndata+(ndata-1)*3;
   xi=xa0(1):(xa0(ndata)-xa0(1))/(nx-1):xa0(ndata);
   
   bxi=interp1(xa0,bxs./b0,xi,'spline');
   byi=interp1(xa0,bys./b0,xi,'spline');
   bzi=interp1(xa0,bzs./b0,xi,'spline');
         
   %--------------------------%
   ny=151;
   %--------------------------%
   disp(' ')
   disp(['nx=',num2str(nx),' ','ny=',num2str(ny)]);
   
   ht=0;
   py=0.1/1;
   mid=round(ny/2)+ht;
   
   x=xi./L0;
   hx=x(2)-x(1);%uniform grids
   hy=py*hx;
   
   disp(['mid=',num2str(mid,'%2.0f'),' ','nYup=',num2str((ny-mid),'%2.0f')]);
   disp(['hx =',num2str(hx,'%2.3f'),' ','hy  =',num2str(hy,'%2.3f')]);
   
   for j=1:ny
      y(j)=(j-mid)*hy;
   end
   disp(['Ymin=',num2str(y(1)*L0,'%2.0f'),' ','Ymax=',num2str(y(ny)*L0,'%2.0f')]);
         
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
      
      close all
      left=0.15; bottom=0.5; width=0.6; height=0.4;
      
      ha=axes('Position',[left bottom width height]);
      
      pcolor(x*L0,y(gg)*L0,Bzup(gg,:)*b0*1e9);
      shading interp;      
      caxis([-60 30]);
      
      hbar=colorbar('vertical');
      pos=get(hbar,'Position');
      set(hbar,'Position',[pos(1)+0.1 pos(2)+0.15 0.01 0.1],'Fontsize',7)
      set(get(hbar,'XLabel'),'String','Bz [nT]','Rotation',0,...
          'HorizontalAlignment','left','Fontsize',8);
      
      hold on
      [cc,hh]=contour(x*L0,y(gg)*L0,Aup(gg,:),[-4:0.11:4],'k');
      set(hh,'linewidth',1.0);
      
      hr1=quiver(xa0,zeros(1,ndata),bxs',bys',0.3,'w');
      set(hr1,'linewidth',1.0);
      
      set(ha,'fontsize',10,'TickDir','out','linewidth',1.0);      
      axis equal
      axis([0 max(x*L0) min(y*L0) max(y*L0)])   
      xlabel('x [km]','fontsize',10)
      ylabel('y [km]','fontsize',10)
    
%=====================================%
% End Original GS Reconstruction Code %
%=====================================%


% --- Outputs from this function are returned to the command line.
function varargout = dve_display_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
axes(handles.ax_bfield);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});
