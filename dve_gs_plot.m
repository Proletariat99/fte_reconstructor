
function[x]=dve_gs_plot(secs,vp,bc,zs,V0,ndata, i1, i2, dzn, Tp)
% Purpose: plots B-field figure using input from dve_angle_loop.
nden=1e6;         %factor for density
nt  =1e6;         %factor for temperature
nv  =1e3;         %factor for velocity
nb = 1e-9;
kb  =1.38*1e-23;  %Boltzmann's constant
miu =4.0*pi*1e-7; %permeability
mp  =1.67*1e-27;  %proton mass

vza=V0*zs';                     % adjusted HT velocity in l ( or z') direction
vzav=vza*zs;                    % adjusted HT velocity in 3D
vhtsv=V0-vzav;                  % 3D velocity in HT frame
vhts=sqrt(vhtsv*vhtsv');		  % magnitude of velocity vector
vhtsn=vhtsv./vhts;              % normalized 3D velocity of plasma in HT frame (in km/s)

xs=-vhtsn;                          % sets xs = to negative 3D velocity of plasma in HT
ys=cross(zs,xs); ys=ys./norm(ys);   % crosses z by x and ensures norm =1
XZ = [xs;ys;zs];                    % makes a new array called XZ
vxs=-V0*xs';                        % x component of V0
% Display x, y, z axes %
disp(['zs =',num2str(zs,'%2.3f  ')]);       % same
disp(['xs =',num2str(xs,'%2.3f  ')]);       % same
disp(['ys =',num2str(ys,'%2.3f  ')]);       % same
disp(['vxs=',num2str(vxs,'%2.1f')]);        % same
disp(' ')

%Project B into reconstruction coord.
bxs=(bc.*nb)*xs';
bys=(bc.*nb)*ys';
bzs=(bc.*nb)*zs';

%Project V into reconstruction coord.
vxc=zeros(1,ndata);
vyc=zeros(1,ndata);
vzc=zeros(1,ndata);
for i=1:ndata
  vxc(i)=(vp(i,:)-V0)*xs';
  vyc(i)=(vp(i,:)-V0)*ys';
  vzc(i)=(vp(i,:)-V0)*zs';
end

%Time to distance
dTim1 = secs(i1:i2)-secs(i1);               % my code (dd, 2011)
xa0=vxs.*dTim1;                             % note: this is NOT deltaX, but sigmax (deltax requires a loop from i1 to i2)

%Calculate vector potential A along y=0
len=length(dTim1);                % # of steps in window = len
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
%-------------------------------------------------------------------------%
%=============================== Block 1 =================================%
%-------------------------------------------------------------------------%
opn=1;
if opn==1
  clf;
  % vmg=sqrt((vxc./vv0).^2 + (vyc./vv0).^2 + (vzc./vv0).^2);    Unused
  % bmg=sqrt((bxs./b0).^2 + (bys./b0).^2 + (bzs./b0).^2);       Unused
  xm=1:ndata;
  
  left=0.1; bottom=0.7; width=0.35; height=0.25; dh=height+0.05; dw=width+0.08;
  
  ha=axes('position',[left bottom width height]);
  plot(xm,vxc./vv0,'k-o',xm,vyc./vv0,'b-o',xm,vzc./vv0,'r-o','Markersize',3); grid on;
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
end             % end for opn ==1

%-------------------------------------------------------------------------%
%=============================== Block 2 =================================%
%-------------------------------------------------------------------------%
opm=1;
if opm==1
  %--interpolation--%
  nx=ndata+(ndata-1)*3;
  xi=xa0(1):(xa0(ndata)-xa0(1))/(nx-1):xa0(ndata);
  
  bxi=interp1(xa0,bxs./b0,xi,'spline');
  byi=interp1(xa0,bys./b0,xi,'spline');
  % bzi=interp1(xa0,bzs./b0,xi,'spline'); % unused here.
  
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
  
  %----------------------------------------------------------------------------%
  %=============================== SubBlock 1 =================================%
  %----------------------------------------------------------------------------%
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
    
    gg=1:ny;                % index for something (1-151)
    Bzup=zeros(ny,nx);      % size of gg by nx (0's 151 by 53)
    for j=1:ny              % Fills in Bzup with data from gsup/gsdn
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
    left=0.15; bottom=0.65; width=0.6; height=0.4;      % parameters for 'Position'
    % left=0.0; bottom=0.0; width=0.6; height=0.4;
    ha=axes('Position',[left bottom width height]);     % creates window in specified area
    B3 = Bzup(gg,:)*b0*1e9;              % defines 3rd axis
    pcolor(x*L0,y(gg)*L0,B3);            % displays 3rd axis (bz) in faceted shading
    shading interp                       % changes shading to interpolated
    %caxis([-60 30]);                     % 3rd axis range (depth/color)
    
    minb3=min(min(B3));
    maxb3=max(max(B3));
    caxis([minb3 maxb3]);
    
    hbar=colorbar('vertical');
    pos=get(hbar,'Position');
    set(hbar,'Position',[pos(1)+0.1 pos(2)+0.1 0.02 0.2],'Fontsize',8)
    set(get(hbar,'XLabel'),'String','Bz [nT]','Rotation',0,...
      'HorizontalAlignment','left','Fontsize',8);
    
    hold on
    %[cc,hh]=contour(x*L0,y(gg)*L0,Aup(gg,:),[-4:0.11:4],'k');  % plots contour lines
    [cc,hh]=contour(x*L0,y(gg)*L0,Aup(gg,:),[-4:0.11:4],'k');
    set(hh,'linewidth',1.0);
    
    hr1=quiver(xa0,zeros(1,ndata),bxs',bys',0.3,'w');  % plots white arrows accross middle.
    set(hr1,'linewidth',1.0);                          % sets linewidth for quiver arrows
    
    set(ha,'fontsize',10,'TickDir','out','linewidth',1.0);
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
    
  end             % end for opm == 2
end               % end for opm == 1
end               % end for function



