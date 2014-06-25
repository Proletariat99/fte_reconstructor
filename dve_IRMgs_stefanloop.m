%GS reconstruction code developed by W.-L. Teh
% Modified by Dave Dyer for THEMIS - 07/2011.
% inputs:   ta = FTE start time -> i1(line #)
%           tb = FTE stop time -> i2(line #)
% outputs:  plot - (B,index);(V,index);(T,index);(N,index) - Line 83ish
%           plot - (x,V) (x,B) (x,A );(x,p_t);(x,dpbz);(A,Bz) (opn == 1)
%           plot - B field lines (opm == 1) & (opm == 2)
%           


function[vht,V0]=dve_IRMgs(doy,ti,tf,sc,source,nev,va_theta,va_phi,theta)
%clear all
%close all
%clc
%------------------%
% HARD CODED STUFF %
%     BEGIN        %
%------------------%
%doy = '2007/06/08';
%ti = '06:30:00';
%tf = '08:30:00';
%sc = 'a';
w = 2              % # of data points on either side of event to add.
nev = 1;           % The event I'm working on (1 is first, not 0) so 1/113, 2/113, etc.
n = 113;           % Total # of FTEs in list (usually just 113)
dhead = 1;         % number of header lines (usually just 1)
%source = 'reduced';% either 'onboard' or 'reduced'.  most are 'onboard'
%va_theta = 0.0;    % theta for variance analysis
%va_phi = 0.0;      % phi for variance analysis
%theta = 0.0;       % theta for local z'
%xs=[0.468420, -0.874707 , -0.124385]
%ys=[0.834420,  0.484265 , -0.263135]
%zs=[0.290400,  0.0194666,  0.956712]

format long;

nden=1e6;         %factor for density
nt  =1e6;         %factor for temperature
nv  =1e3;         %factor for velocity
nb  =1e-9;        %factor for magnetic field
kb  =1.38*1e-23;  %Boltzmann's constant
kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
miu =4.0*pi*1e-7; %permeability
mp  =1.67*1e-27;  %proton mass

%------------------%
% HARD CODED STUFF %
%       END        %
%------------------%
%------------%
% TIME STUFF %
% START      %
%------------%
yr = sscanf(doy,'%4c%*6c',2);
month = sscanf(doy, '%*5c%2c%*3c',2);
day = sscanf(doy, '%*8c%2c',2);
HHi = sscanf(ti,'%2c%0.2*4c');
MMi = sscanf(ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(ti, '%*6c%2c0', 2);
HHf = sscanf(tf,'%2c%0.2*4c');
MMf = sscanf(tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(tf, '%*6c%2c0', 2);

ff=['/FTE/2011/findings/themis_fte_iaga_2011_clean_final_timeorder_recons.txt'];
fid = fopen(ff);
f = textscan(fid, '%2c %s %s %s %f %s %s %f %f %f', n);
fclose(fid);

fn = ['th',sc,'_',num2str(strcat(yr,month,day)),'_',num2str(strcat(HHi,MMi)),'_',num2str(strcat(HHf,MMf)),'_GSM_v4_',source,'.dat'];
a = importdata(fn,' ',dhead);
sec = a.data(:,1);
BX=a.data(:,2); BY=a.data(:,3); BZ=a.data(:,4); 	% Assigns Data columns to Bx, By, Bz arrays (nT(GSE))
VX=a.data(:,6); VY=a.data(:,7); VZ=a.data(:,8); 	% assigns data columns to Vx, Vy, Vz arrays (km/s(GSE))
Nx=a.data(:,10); 					% assigns column 5 (plasma density) to Nx (density, in units cm^-3)		
Tx=a.data(:,11); 					% assigns column 6 (plasma temp in K) to Tx
idata=1:length(BX);				% idata = one row by (# of rows in Bdata)columns

nev = int16(nev);
evn = nev;
sc = f{1}(nev);  % This shit is broken (f{2}(1) = f{2}(nev) when nev = 1
ta = f{2}(nev);
tb = f{4}(nev);
sdoy = strcat(num2str(yr),'-',num2str(month),'-',num2str(day));
tstart = strcat(sdoy,'/',ti);
tstop = strcat(sdoy,'/',tf);
ts1 = datestr(tstart, 'mmmm/dd/yyyy HH:MM:SS:FFF');
ts2 = datestr(tstop, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vts1 = datevec(ts1, 'mmmm/dd/yyyy HH:MM:SS:FFF');
vts2 = datevec(ts2, 'mmmm/dd/yyyy HH:MM:SS:FFF');
%dts = datenum(vts2-vts1)*3600*24
fvts1 = datevec(ta, 'yyyy-mm-dd/HH:MM:SS');
fvts2 = datevec(tb, 'yyyy-mm-dd/HH:MM:SS');
dn1 = datenum(vts1);
dn2 = datenum(vts2);
fdn1 = datenum(fvts1);
fdn2 = datenum(fvts2);

t1f1 = etime(fvts1,vts1);   % returns number of seconds from start of window to start of FTE
t1f2 = etime(fvts2,vts1);   % elapsed time from start of window to stop of FTE (s)
disp(['time window is from ',num2str(datestr(fvts1)),' to ', num2str(datestr(fvts2))])
disp(' ')
disp(['This should be equivalent to line # ', num2str(evn), ' in themis_fte_iaga_2011_clean_final_timeorder_recons.txt'])

tstep = a.data(:,1);
j = tstep >= t1f1 & tstep < t1f2;
k = find(j);
%fwin = a.data(tstep >= t1f1 & tstep < t1f2,1);
ns1 = k(1)-w;                % adds two data points to either side
ns2 = k(end) + w;            % of FTE
ind = [ns1:ns2];             
ind = ind';
%------------%
% TIME STUFF %
% STOP       %
%------------%



%---------------------------------------------%
% Convert t1 % t2 into line # in the GSM file %
%---------------------------------------------%
i1 = ns1;
i2 = ns2;

%Data Plot
subplot(3,1,1);                                                  % creates a subplot
plot(idata,BX,'k',idata,BY,'b',idata,BZ,'r','linewidth',1.0);	% plots idata vs. Bx, By, Bz
subplot(3,1,2);                                                  % creates a subplot
plot(idata,VX,'k',idata,VY,'b',idata,VZ,'r','linewidth',1.0);	% plots idata vs. Vx, Vy, Vz
subplot(3,1,3);                                                  % creates a subplot
plot(idata,Nx,'k',idata,Tx,'r','linewidth',1.0);                % plots idata (x) vs. density & temp


%----------------------%
%HT and MVAB analysis  %
%----------------------%
bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];				% New B matrix x3
vp=[VX(i1:i2) VY(i1:i2) VZ(i1:i2)];				% New V matrix x3
Nc=Nx(i1:i2);                                   % New Density matrix x3

[EHT,EC,HTR,HTslope,vht]=HTcoef(bc,vp,length(bc(:,1)));		% Sets EHT, EC, HTR, Htslope & vht from another program. 

disp(['VHT= ',num2str(vht,'%2.1f  ')]);			% displays Hoffmann-Teller Velocity
disp(' ')
disp(['[cc_HT,slope]=',num2str([HTR HTslope],'%2.2f ')]);	% displays HTR & HTslope
disp(' ');                                      % space
clear('bc','vp','Nc');                          % clears bc, vp and nc variables

bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];				% re-assigns bc for some reason.
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(bc,length(bc(:,1)));	% uses MVAB.m to define variables:
								% input:
								  % bc - bfield 3D array
								  % ncases - number of cases in bc(i1:i2)
								% output:
								  % Kbmax - Max variance?
								  % Kbint - Intermediate variance?
								  % Kbmin - Min variance?
								  % qmva - quality ratio?
								  % d1 - lambda1
								  % d2 - lambda2
								  % d3 - lambda3

Kbint=-Kbint;                                   % sets intermediate to inverse
disp(['d=',num2str([d1 d2 d3]/d3,'%2.1f  ')]);	% no idea.  wtf is d3?
disp(['L= ',num2str(Kbmax,'%2.3f  ')]);			% L coordinate
disp(['M= ',num2str(Kbint,'%2.3f  ')]);			% M coordinate
disp(['N= ',num2str(Kbmin,'%2.3f  ')]);			% Normal coordinate (to what?)
disp(' ');
%------------------------------------
clear('bc');

i1=ns1;                           % start record for reconstruction
i2=ns2;                          % stop record for reconstruction

iix=i1:i2;                      % sets up array of i1 to i2 (columns)		
ndata=length(iix);              % counts iix (i2 - i1)

%Initial Data for Reconstruction
bc=[BX(iix) BY(iix) BZ(iix)];   % bc reconstruction by matrix
vp=[VX(iix) VY(iix) VZ(iix)];	% vp reconstruction by matrix

dzn=Nx(iix);                    % density reconstruction by matrix
%Tp =Tx(iix);                    % temperature reconstruction by matrix
Tp = Tx(iix)/kb2;

opm=1;                          % gate for reconstruction: switch for later 'if' statement.

phi1=va_phi*pi/180;                 % hard coded phi offset (horiz. offset from normal?)
CM=[1 0 0; 0 cos(phi1) -sin(phi1) ; 0 sin(phi1) cos(phi1)];		% 3x3 identity matrix (if phi = 0)
XZ=[Kbmax;Kbint;Kbmin];			% 3x3 matrix with MVAB output as rows

disp(['XZ 0= '])
disp(num2str(XZ))
disp(' ')

XZ=CM*XZ;                       % matrix multiplication with phi

disp(['XZ 1= '])
disp(num2str(XZ))
disp(' ')

theta1=va_theta*pi/180;              % hard coded theta (40?) offset (vert. offset from invariant axis)
CM=[cos(theta1) -sin(theta1) 0 ; sin(theta1) cos(theta1) 0 ; 0 0 1];% 3x3 identify matrix (if theta = 0)
XZ=[XZ(1,:); XZ(2,:) ; XZ(3,:)];% adjusts XZ for phi (XZ = XZ if CM = identity matrix)
XZ=CM*XZ;                       % Adjusts XZ for theta (XZ = XZ if CM = identity matrix)

disp(['XZ 2= '])
disp(num2str(XZ))
disp(' ')

zs=XZ(2,:);                     % invariant axis (l or z')

V0 = [-60.0000, 94.3000, -6.10000];    % SHEATH VALUE
%V0 = [-55.1000, 90.5913, 14.8534];  % PAPER VALUE
zs=[-0.516, -0.301, -0.802];
%V0=vht;                         % V0 - reference frame velocity (these can be different, actually)

%%%%%%%% HERE STARTS STEFAN LOOPS FOR zs %%%%%%%%
pol=5.0;                 % in deg
azi=20.0;                % in deg (make variable in the input)

uinput='';               %formerly known as 'ans'
nxorig=zs(1);
nyorig=zs(2);
nzorig=zs(3);

theta0=acosd(nzorig);
  phi0=atand(nyorig/nxorig);  
  
if (nxorig < 0.0) && (nyorig >= 0.0);
    phi0=-phi0+pi;
end
if (nxorig < 0.0) && (nyorig < 0.0); 
    phi0= phi0+pi;
end
if (nxorig >= 0.0) && (nyorig < 0.0); 
    phi0=-phi0+pi*2.;
end

disp(['phi0=',num2str(phi0)])
disp(['theta0=',num2str(theta0)])
%nxtest=sin(theta0)*cos(phi0)
%nytest=sin(theta0)*sin(phi0)
%nztest=cos(theta0)
%print,nxorig,nxtest
%print,nyorig,nytest
%print,nzorig,nztest

%theta0=theta0*180./pi;
%  phi0=phi0*180./pi;
%print,'Theta0 (deg): ',theta0
%print,'Phi0 (deg): ',phi0

num1=(round(90./pol)*2);
num2=round(360./azi);
ntheta=zeros(num1,1);
for i=1:2:(num1-1); 
    ntheta(i,1) = (theta0-(i/2)*pol);
end;
for i=2:2:num1;
    ntheta(i,1)=theta0+(i+1)/2*pol;
end;
    
nphi=zeros(num2,1);
for i=1:1:num2
    nphi(i,1)=phi0+i*azi
end
for j=1:num1 
  theta=ntheta(j)
  for k=1:num2
       phi=nphi(k)
       if theta > 360.0 
           theta=theta-360.0
       end
       if theta < -360.0 
           theta=theta+360.0
       end
       if phi > 360.0 
           phi=phi-360.0
       end
       nx=sin(theta*pi/180.)*cos(phi*pi/180.)
       ny=sin(theta*!dpi/180.)*sin(phi*!dpi/180.)
       nz=cos(theta*!dpi/180.)
       nitr=[nx,ny,nz]
       nmag=sqrt(total(nitr*nitr))
       print,'Theta new and old (deg): ',theta,theta0
       print,'Phi new and old (deg): ',phi,phi0
       print,'GSM unit vector: ',nitr(0),nitr(1),nitr(2),nmag

       zs=nitr
%%%%%%%% NOW BACK TO ORIGINAL MATLAB ROUTINES FROM STEFAN IDL LOOP %%%%%%%%

vza=V0*zs';                     % adjusted HT velocity in l ( or z') direction
%vza = vht*zs';
vzav=vza*zs;                    % adjusted HT velocity in 3D
vhtsv=V0-vzav;                  % 3D velocity in HT frame
%vhtsv=vht-vzav;
vhts=sqrt(vhtsv*vhtsv');		% no idea - magnitude of velocity vector, i guess
vhtsn=vhtsv./vhts;              % normalized 3D velocity of plasma in HT frame (in km/s)

xs=-vhtsn;                      % sets xs = to negative 3D velocity of plasma in HT
ys=cross(zs,xs); ys=ys./norm(ys);   % finds y' by crossing x' & z' and normalizing to 1

% Display Min Variance Results
disp('XZ Matrix is:')
disp(XZ)
disp(['ZS is ',num2str(zs)])
disp(' ')
disp(['Kbmin = ',num2str(Kbmin)])
disp(['Kbint = ',num2str(Kbint)])
disp(['Kbmax = ',num2str(Kbmax)])
disp(' ')

theta2=theta*pi/180;
CM=[cos(theta2) -sin(theta2) 0 ; sin(theta2) cos(theta2) 0 ; 0 0 1];
XZ=[xs;ys;zs];
XZ=CM*XZ;

xs=XZ(1,:); xs=xs/norm(xs);     % Note: norm(xs) = magnitude of xs
ys=XZ(2,:); ys=ys/norm(ys);     % adjusts again for theta (above)

vxs=-V0*xs';                    % re-adjust V0 by theta

disp(['ndata=',num2str(ndata,'%2.0f')]);    % displays results
disp(['zs =',num2str(zs,'%2.3f  ')]);       % same
disp(['xs =',num2str(xs,'%2.3f  ')]);       % same
disp(['ys =',num2str(ys,'%2.3f  ')]);       % same
disp(['vxs=',num2str(vxs,'%2.1f')]);        % same
disp(' ')

bx=bc(:,1);	by=bc(:,2);	bz=bc(:,3);         % bx, by, bz
vx=vp(:,1);	vy=vp(:,2);	vz=vp(:,3);         % vx, vy, vz

%Project B into reconstruction coord.
bxs=(bc.*nb)*xs';                           
bys=(bc.*nb)*ys';
bzs=(bc.*nb)*zs';

%Project V into reconstruction coord.
for i=1:ndata
   vxc(i)=(vp(i,:)-V0)*xs';
   vyc(i)=(vp(i,:)-V0)*ys';
   vzc(i)=(vp(i,:)-V0)*zs';   
end

dTim1 = sec(i1:i2)-sec(i1);               % my code
%Time to distance
xa0=vxs.*dTim1;                             % note: this is NOT deltaX, but sigmax (deltax requires a loop from i1 to i2)
disp(['first 5 time values are '])
disp(num2str(dTim1(1:5)))
disp(' ')

%Calculate vector potential A along y=0
ndA(1)=0;
len=length(dTim1);
for i=2:len
   ndA(i)=-(bys(i)+bys(i-1))*(dTim1(i)-dTim1(i-1))*vxs*nv*0.5;
end
A1=cumsum(ndA);                             % Not sure if we still need this line.
A1=A1';

bb=sqrt(bxs.^2+bys.^2+bzs.^2);
%pp=dzn.*Tp*kb*nden*nt;                 % old pp factor
pp=dzn.*Tp*kb*nden;

bmax=max(bb);
nmax=max(dzn)*nden;

%Normalized factors
b0=bmax;
p0=bmax*bmax/miu;
n0=nmax;
T0=p0/(kb*n0);
v0=sqrt(b0^2/(miu*n0*mp)*1e-6);
A0=max(abs(A1));
L0=(A0/bmax)*1e-3;
pbz=pp./p0+((bzs./b0).^2)/2;
An=A1./A0;

opn=1;
if opn==1
   clf;
   vmg=sqrt((vxc./v0).^2 + (vyc./v0).^2 + (vzc./v0).^2);
   bmg=sqrt((bxs./b0).^2 + (bys./b0).^2 + (bzs./b0).^2);
   xm=1:ndata;
   
   left=0.1; bottom=0.7; width=0.35; height=0.25; dh=height+0.05; dw=width+0.08;
   
   ha=axes('position',[left bottom width height]);
   plot(xm,vxc./v0,'k-o',xm,vyc./v0,'b-o',xm,vzc./v0,'r-o','Markersize',3); grid on; 
   ylabel('V','rotation',0);
   set(ha,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
   %--
   ha=axes('position',[left bottom-dh width height]);
   plot(xm,bxs./b0,'k-o',xm,bys./b0,'b-o',xm,bzs./b0,'r-o','Markersize',3); grid on; 
   ylabel('B','rotation',0);
   set(ha,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
   %--
   ha=axes('position',[left bottom-dh*2 width height]);
   plot(An,'k-o','Markersize',3);  grid on; 
   ylabel('A','rotation',0);
   set(ha,'ylim',[-1 1],'ytick',[-1:0.25:1],'fontsize',8);
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
   
   ha=axes('position',[left+dw bottom width height]);
   plot(An,pbz,'ko','Markersize',5); grid on;    
   hold on;
   plot(ax1,yys1,'color',[0.5 0.5 0.5],'linewidth',3);    
   set(ha,'ylim',[min(pbz)-0.1 max(pbz)+0.1],'ytick',[0:0.1:1],'fontsize',8);
   ylabel('p_t','rotation',0);
   xlabel('A')
   %--
   ha=axes('position',[left+dw bottom-dh width height]);
   plot(ax1,polyval(df1,ax1),'k'); grid on;
   ylabel('Dpbz');
   set(ha,'fontsize',8);
   %--
   ha=axes('position',[left+dw bottom-dh*2 width height]);
   plot(An,bzs/b0,'ko','Markersize',3); grid on; 
   hold on;
   plot(ax1,yyz1,'b','linewidth',3);
   set(ha,'ylim',[min(bzs/b0)-0.1 max(bzs/b0)+0.1],'ytick',[-1:0.25:1],'fontsize',8);
   ylabel('Bz');
   xlabel('A')
end


if opm==1
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
   
   opn=2;      
   if opn==2
      
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
      %axis([0 max(x*L0) -500 500])`
      xlabel('x [km]','fontsize',10)
      ylabel('y [km]','fontsize',10)
      
   end

%%%%%%%% RE-ENTER STEFAN ADDED LOOPS %%%%%%%%
       read,'Enter (0/1):',ans
       if ans ne '' then goto, stop
       i=i+1L
  endfor
endfor
stop:
%%%%%%%% EXIT STEFAN ADDED LOOPS %%%%%%%%

end

