%GS reconstruction code developed by W.-L. Teh
%modified for use with THEMIS by D. Dyer - June, 2011

clear all
close all
clc

format long;

nden=1e6;         	%factor for density
nt  =1 			%new factor for temperature (1)
%1e6;        		% old factor for temperature
nv  =1e3;	  	%factor for velocity
nb  =1e-9;        	%factor for magnetic field

kb = 8.61739*1e-5;	%Boltzmann's Constant (in eV/K)
%kb  =1.38*1e-23;  	%Boltzmann's constant (in mK)
miu =4.0*pi*1e-7; 	%permeability
mp  =1.67*1e-27;  	%proton mass

%Magnetic Field Input
%str1='Bdata.dat';
str1 = 'thc_20080415_0600_1000_GSM_v4_onboard.dat'
fid1=fopen(str1,'rt')
fgets(fid1)
Data1=fscanf(fid1,'%f %f%f%f',[4,inf]);
fclose(fid1);
Data1=Data1';

secb=Data1(:,1); 				%sec
BX=Data1(:,2); BY=Data1(:,3); BZ=Data1(:,4); 	%nT(GSE)

clear('Data1');

%Plasma Data Input
str1='Pdata.dat'

fid1=fopen(str1,'rt');
fgets(fid1);
Data1=fscanf(fid1,'%f %f%f%f %f %f',[6,inf]);
fclose(fid1);
Data1=Data1';

secp=Data1(:,1); %sec
VX=Data1(:,2); VY=Data1(:,3); VZ=Data1(:,4); %km/s(GSE)
Nx=Data1(:,5); %/cc
Tx=Data1(:,6); %MK

clear('Data1');

idata=1:length(BX);

%Data Plot
subplot(3,1,1)
plot(idata,BX,'k',idata,BY,'b',idata,BZ,'r','linewidth',1.0);
subplot(3,1,2)
plot(idata,VX,'k',idata,VY,'b',idata,VZ,'r','linewidth',1.0);
subplot(3,1,3)
plot(idata,Nx,'k',idata,Tx,'r','linewidth',1.0);


%HT and MVAB analysis
i1=1;
i2=27;

bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];
vp=[VX(i1:i2) VY(i1:i2) VZ(i1:i2)];
Nc=Nx(i1:i2);

[EHT,EC,HTR,HTslope,vht]=HTcoef(bc,vp,length(bc(:,1)));

disp(['VHT= ',num2str(vht,'%2.1f  ')]);
disp(['[cc_HT,slope]=',num2str([HTR HTslope],'%2.2f ')]);
disp(' ');
clear('bc','vp','Nc');

bc=[BX(i1:i2) BY(i1:i2) BZ(i1:i2)];
[Kbmax,Kbint,Kbmin,qmva,d1,d2,d3]=MVAB(bc,length(bc(:,1)));

Kbint=-Kbint;
disp(['d=',num2str([d1 d2 d3]/d3,'%2.1f  ')]);
disp(['L= ',num2str(Kbmax,'%2.3f  ')]);
disp(['M= ',num2str(Kbint,'%2.3f  ')]);
disp(['N= ',num2str(Kbmin,'%2.3f  ')]);
disp(' ');
%------------------------------------
clear('bc');


i1=1;
i2=30;

iix=i1:i2;
ndata=length(iix);

%Initial Data for Reconstruction
bc=[BX(iix) BY(iix) BZ(iix)];
vp=[VX(iix) VY(iix) VZ(iix)];

dzn=Nx(iix);
Tp =Tx(iix);

opm=1; %gate for reconstruction


phi=0.0*pi/180;
CM=[1 0 0; 0 cos(phi) -sin(phi) ; 0 sin(phi) cos(phi)];
XZ=[Kbmax;Kbint;Kbmin];
XZ=CM*XZ;

theta=40.0*pi/180;
CM=[cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1];
XZ=[XZ(1,:); XZ(2,:) ; XZ(3,:)];
XZ=CM*XZ;

zs=XZ(2,:) %invariant axis
V0=vht %reference frame velocity

vza=V0*zs';
vzav=vza*zs;
vhtsv=V0-vzav;
vhts=sqrt(vhtsv*vhtsv');
vhtsn=vhtsv./vhts;

xs=-vhtsn;
ys=cross(zs,xs); ys=ys./norm(ys);

theta=0.0*pi/180;
CM=[cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1];
XZ=[xs;ys;zs];
XZ=CM*XZ;

xs=XZ(1,:); xs=xs/norm(xs);
ys=XZ(2,:); ys=ys/norm(ys);

vxs=-V0*xs';

disp(['ndata=',num2str(ndata,'%2.0f')]);
disp(['zs =',num2str(zs,'%2.3f  ')]);
disp(['xs =',num2str(xs,'%2.3f  ')]);
disp(['ys =',num2str(ys,'%2.3f  ')]);
disp(['vxs=',num2str(vxs,'%2.1f')]);


bx=bc(:,1);	by=bc(:,2);	bz=bc(:,3);
vx=vp(:,1);	vy=vp(:,2);	vz=vp(:,3);

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

[cx,cy]=size(secb);
for i=2:cx
   dTime(1,i)=secb(1);
end
dTime(1,1)=0;
dTim1=cumsum(dTime)

%Time to distance
xa0=vxs.*dTim1;

%Calculate vector potential A along y=0
ndA(1)=0;
len=length(dTim1);
for i=2:len
   ndA(i)=-(bys(i)+bys(i-1))*(dTim1(i)-dTim1(i-1))*vxs*nv*0.5;
end
A1=cumsum(ndA);
A1=A1';

bb=sqrt(bxs.^2+bys.^2+bzs.^2);
pp=dzn.*Tp*kb*nden*nt;

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
      xlabel('x [km]','fontsize',10)
      ylabel('y [km]','fontsize',10)
      
   end

end

