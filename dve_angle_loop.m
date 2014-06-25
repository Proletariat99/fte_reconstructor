% Program name: dve_gs_plot.m
% Purpose: To manually iterate through 360 degrees
% Inputs: doy - date in format:'2011-07-18'
%         ti - start time of KH window in format:'2011-07-18/09:52:16'
%         tf - stop time of KH window in format:'2011-07-18/09:53:32'
%         sc - THEMIS spacecraft designation in format:'a', 'b', or 'c'
%         src - THEMIS Electrostatic Analyzer (ESA) source in format:'reduced' or
%         'onboard'
%         nev - line number of flux transfer event (FTE) in text file in
%         format: 01 (no quotes)
%         n - total number of FTEs in text file in format: 113 (no quotes)
%         w - # of data points to add on either side of FTE window in format: 02 (no
%         quotes)
%         FROM FILES:
%           event start
% Outputs:  1) Plots of A vs. various things (optional?
%           2) B field line plots
%           3) Switch for "accept" or "reject and move to next theta/phi"
%           4)
% Author:   Dave Dyer, July, 2011
% Changelog:    2011-07-18: Initial Program Creation.

function[zs]=dve_angle_loop(doy, ti, tf, sc, src, nev, n, w, pol, azi)
% Purpose: To manually iterate through 360 degrees

%-----------%
% Constants %
%-----------%
nden=1e6;         %factor for density
nt  =1e6;         %factor for temperature
nv  =1e3;         %factor for velocity
nb  =1e-9;        %factor for magnetic field
kb  =1.38*1e-23;  %Boltzmann's constant
kb2 = 8.61739*1e-5;   % Boltzmann's constant (eV/K);
miu =4.0*pi*1e-7; %permeability
mp  =1.67*1e-27;  %proton mass
dhead = 1;        % # of header lines in time_step file

%-------------%
% Input Check %
%-------------%
disp(['doy=',doy])
disp(['ti=',ti])
disp(['tf=',tf])
disp(['sc=',sc])
disp(['source=',src])
disp(['nev=',num2str(nev)])
disp(['n=',num2str(n)])
disp(['w=',num2str(w)])

%-----------------------------------%
% Read Magnetosheath File for Input %
%-----------------------------------%
fnsh=['/FTE/2011/findings/themis_fte_iaga_2011_msheath.txt'];
mshfid = fopen(fnsh);
msh = textscan(mshfid, '%s %s %s %s %f %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', n);
fclose(mshfid);

%-----------------------------------%
% An error check for a small mistake%
% that can cause a lot of problems  %
%-----------------------------------%
nevtest=ischar(nev);
ntest=ischar(n);
wtest=ischar(w);
tests = [nevtest,ntest,wtest];
errtest = any(tests);
switch(errtest)
  
  case 1
    disp(' ')
    disp('============================================')
    disp('============================================')
    disp('CHANGE nev,n or w TO A FLOATING POINT VALUE!')
    disp('============================================')
    disp('============================================')
    disp(' ')
    disp(['nev=',class(nev)])
    disp(['n=',class(n)])
    disp(['w=',class(w)])
    return
    %-----------------%
    % end error check %
    %-----------------%
  otherwise
end             % end errortest
%-------------------------------%
% set xs to V-magnetosheath (V0)%
%-------------------------------%
msh_vx = msh{13}(nev);
msh_vy = msh{14}(nev);
msh_vz = msh{15}(nev);
msh_v = [msh_vx,msh_vy,msh_vz];

disp(['msh_vx=',num2str(msh_vx)])
disp(['msh_vy=',num2str(msh_vy)])
disp(['msh_vz=',num2str(msh_vz)])

% V0=msh_v/norm(msh_v);             Don't normalize V0 for distance calculations.
V0=msh_v
nV0=V0/norm(V0)
xs = -nV0;

disp(['V0 = ',num2str(V0)])
disp(['xs = ',num2str(xs)])

%------------%
% time stuff %
% start      %
%------------%
tcf = 24*3600.;  %time conversion factor
yr = sscanf(doy,'%4c%*6c',2);
month = sscanf(doy, '%*5c%2c%*3c',2);
day = sscanf(doy, '%*8c%2c',2);
HHi = sscanf(ti,'%2c%0.2*4c');
MMi = sscanf(ti, '%*3c%2c0%*2c', 2);
SSi = sscanf(ti, '%*6c%2c0', 2);
HHf = sscanf(tf,'%2c%0.2*4c');
MMf = sscanf(tf, '%*3c%2c0%*2c', 2);
SSf = sscanf(tf, '%*6c%2c0', 2);

%----------------------%
% import FTE file to f %
%----------------------%
ff='/FTE/2011/findings/themis_fte_iaga_2011_clean_final_timeorder_recons.txt';
fid = fopen(ff);
f = textscan(fid, '%2c %s %s %s %f %s %s %f %f %f', n);
fclose(fid);

%---------------------------------%
% import time step data file to a %
%---------------------------------%
fn = ['th',sc,'_',num2str(strcat(yr,month,day)),'_',num2str(strcat(HHi,MMi)),'_',num2str(strcat(HHf,MMf)),'_GSM_v4_',src,'.dat'];
a = importdata(fn,' ',dhead);
secs = a.data(:,1);
BX=a.data(:,2); BY=a.data(:,3); BZ=a.data(:,4); 	% Assigns Bx, By, Bz (from file - GSM)
VX=a.data(:,6); VY=a.data(:,7); VZ=a.data(:,8); 	% Assigns Vx, Vy, Vz (from file - (km/s(GSE))
Nx=a.data(:,10); 					% assigns plasma N density from file (cm^-3)
Tx=a.data(:,11); 					% assigns temp from file (K)
idata=1:length(BX);                 % index from 1 to (# of BX)

sc = f{1}(nev);
sdoy = strcat(num2str(yr),'-',num2str(month),'-',num2str(day));
tstart = strcat(sdoy,'/',ti);
tstop = strcat(sdoy,'/',tf);
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
disp(['Data Window Start = ',num2str(ti)])
disp([' Data Window Stop = ',num2str(tf)])
disp(' ')
%---------------------------------%
% 1 & 2 consider stop/start times %
% for FTE from the file or array  %
%---------------------------------%
t1 = f{2}(nev);
t2 = f{4}(nev);
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
tstep = a.data(:,1);
j = tstep >= t1f1 & tstep < t1f2;
k = find(j);
s1 = (k(1)-w);
s2 = (k(end)+w);
ftdnd = secs(s1)/(3600*24);       % FTE Window Time Down (in date number)
ftdnu = secs(s2)/(3600*24);       % FTE Window Time Up (in date number)
ftdndadj = dni+ftdnd;
ftdnuadj = dni+ftdnu;
fdna = ftdndadj;                  % move up/down
fdnb = ftdnuadj;                  % talk about an unwieldy variable...
fvta = datevec(fdna);
fvtb = datevec(fdnb);
ns1 = k(1)-w;                % adds two data points to either side
ns2 = k(end) + w;            % of FTE
ind = [ns1:ns2];
ind = ind';
ta = datestr(fvta);
tb = datestr(fvtb);
disp('===========================Additional Data Points=========================')
disp(['w=',num2str(w)])
disp(['Start Line (file): ',num2str(k(1))])
disp(['  End Line (file): ' ,num2str(k(end))])
disp(['Start Line (manually adjusted): ',num2str(ns1)])
disp(['  End Line (manaully adjusted): ' ,num2str(ns2)])
disp(' ')
disp(['Start (file): ',datestr(fdn1)])
disp(['  End (file): ',datestr(fdn2)])
disp(['              (', num2str((fdn2-fdn1)*tcf),' total seconds.)'])
disp(' ')
disp(['Start (manually adjusted): ',datestr(fvta)])
disp(['  End (manually adjusted): ',datestr(fvtb)])
disp(['                           (', num2str((fdnb-fdna)*tcf),' total seconds.)'])
disp(' ')
disp('================================End Data Points============================')

%-------------%
% time stuff  %
% stop        %
%-------------%

%------------------%
% get zs from MVAB %
%------------------%
i1 = ns1;
i2 = ns2;
bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];				% Bfield from file (GSM)
vp=[VX(i1:i2) VY(i1:i2) VZ(i1:i2)];				% V from file (GSM)
Nc=Nx(i1:i2);                                   % Density for window ta to tb

[EHT,EC,HTR,HTslope,vht]=HTcoef(bc,vp,length(bc(:,1)));

clear('bc','vp','Nc');                          % clears bc, vp and nc variables


bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];				% re-assigns B for some reason, should be the same as above.

[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(bc,length(bc(:,1)));	% uses MVAB.m to define variables:
Kbint=-Kbint;                                   % z = -z
disp(['d=',num2str([d1 d2 d3]/d3,'%2.1f  ')]);	
disp(['L(X)= ',num2str(Kbmax,'%2.3f  ')]);			% X
disp(['M(Z)= ',num2str(Kbint,'%2.3f  ')]);			% Z
disp(['N(Y)= ',num2str(Kbmin,'%2.3f  ')]);			% Y
disp(' ');

zs = Kbint;
ys = cross(zs,xs);

disp('================Local Coordinate System (GSM)==================')
disp(['xs=',num2str(xs)])
disp(['ys=',num2str(ys)])
disp(['zs=',num2str(zs)])
disp('================Local Coordinate System (GSM)==================')
disp(' ')
% Display Min Variance Results
disp(['Kbmin = ',num2str(Kbmin)])
disp(['Kbint = ',num2str(Kbint)])
disp(['Kbmax = ',num2str(Kbmax)])
disp(' ')
iix=i1:i2;                      % sets up array of i1 to i2 (columns)
disp('=====================MVAB Results (above)=======================')
ndata=length(iix);              % counts iix (i2 - i1)
%Initial Data for Reconstruction
bc=[BX(iix) BY(iix) BZ(iix)];   % bc reconstruction by matrix
vp=[VX(iix) VY(iix) VZ(iix)];	% vp reconstruction by matrix

dzn=Nx(iix);                    % density reconstruction by matrix
Tp = Tx(iix)/kb2;

%-----------------------------------------------------------------------------------%
%===========================Begin Stefan Loop Area==================================%
%-----------------------------------------------------------------------------------%
nxorig=zs(1);
nyorig=zs(2);
nzorig=zs(3);
theta0=acosd(nzorig);
phi0=atand(nyorig/nxorig);
%-------------------------------------%
% set phi0 based on nxorig and nyorig %
%-------------------------------------%
if (nxorig < 0.0) && (nyorig >= 0.0);
  phi0=-phi0+180.;
end
if (nxorig < 0.0) && (nyorig < 0.0);
  phi0= phi0+180.;
end
if (nxorig >= 0.0) && (nyorig < 0.0);
  phi0=-phi0+360.;
end
disp(['phi0=',num2str(phi0)]);
disp(['theta0=',num2str(theta0)]);

num1=round(180/pol)*2;
% %-----------------------------%
% % Error Check for polar angle %
% %-----------------------------%
% if (180/pol) < 1
%   disp('GO BACK AND CHOOSE A SMALLER POLAR ANGLE. (<90)')
%   return;
% end
% %-----------------%
% % END error check %
% %-----------------%
ntheta=zeros((num1),1);    % initialize ntheta (with zeros)

for ii=1:2:(num1-1);         % populates ntheta (odd)
  ntheta(ii,1) = (theta0-(ii-1)*pol/2);
end;                        
for jj=2:2:num1;             % populates ntheta (even)
  ntheta(jj,1)=((theta0+(jj-1)*pol/2)+pol/2);
end;
%disp(['ntheta is', num2str(ntheta(1:10)')]);

num2=round(360/azi);
nphi=zeros(num2,1); % pre-allocate array.
%nphi(1,1)=phi0;     % first value is phi0 itself.
for kk=1:1:num2      % populates nphi (all)
  nphi(kk,1)=phi0+(kk-1)*azi;
end;   
disp(['nphi is ',num2str(nphi')]);

%==================%
% BIG LOOP (theta) %
%==================%
for j=1:num1
  theta=ntheta(j);
  %================%
  % BIG LOOP (phi) %
  %================%
  for k=1:num2
    phi=nphi(k);
    %----------------------%
    % Check for Big Angles %
    %----------------------%
    if theta > 360.0
      theta=theta-360.0;
    end
    if theta < -360.0
      theta=theta+360.0;
    end
    
    if phi > 360.0
      phi=phi-360.0;
    end
    %--------------------------%
    % End Check for Big Angles %
    %--------------------------%
    nx=sind(theta)*cosd(phi);
    ny=sind(theta)*sind(phi);
    nz=cosd(theta);
    nitr=[nx,ny,nz];
    disp('----------angles---------')
    disp(['Theta old= ',num2str(theta), ''])
    disp(['Theta new= ',num2str(theta0), '?'])
    disp(['Phi old = ',num2str(phi), '?'])
    disp(['Phi new = ',num2str(phi0), '?'])
    disp('----------angles---------')
    disp(' ')
    
    zs=nitr;
    
    [dTim1]=dve_gs_plot(secs,vp,bc,zs,V0,ndata, i1, i2, dzn, Tp);
    % dve_angledata contains confirmation dialog and labels.
    [zs]=dve_angledata(pol,azi,theta,phi,phi0, theta0,zs,V0,j, k, d1 ,d2, d3); 
  end;
  %==============================%
  % END INNER LOOP (phi - k end) %
  %==============================%  
end;
%================================%
% END OUTER LOOP (theta - j end) %
%================================%



