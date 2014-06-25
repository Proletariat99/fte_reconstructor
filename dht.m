%Dehoffmann Teller Frame Velocity
%B in nT, V in km/s
function [vht]=dht(bc,vp,ncase)
 va0=100./((4.*pi*1.67)^0.5);  
 

 bx=bc(:,1);
 by=bc(:,2);
 bz=bc(:,3);
  
 vx=vp(:,1);
 vy=vp(:,2);
 vz=vp(:,3);
 


        b2=(bx.^2.0+by.^2.0+bz.^2.0);
        bv=(bx.*vx+by.*vy+bz.*vz);
 
%  **   calculate (K.v)o -- a vector
        kvx=0; kvy=0; kvz=0;
        for j=1:ncase
        kvx=kvx+b2(j)*vx(j)-bx(j)*bv(j);
        kvy=kvy+b2(j)*vy(j)-by(j)*bv(j);
        kvz=kvz+b2(j)*vz(j)-bz(j)*bv(j);
        end
        kvxs=kvx; kvys=kvy; kvzs=kvz;
        kvx=kvx/ncase;
        kvy=kvy/ncase;
        kvz=kvz/ncase;
 
%  **   calculate Ko -- a matrix
        k11=0; k12=0; k13=0;
        k21=0; k22=0; k23=0;
        k31=0; k32=0; k33=0;
        for j=1:ncase
        k11=k11+b2(j)-bx(j)*bx(j);
        k12=k12-bx(j)*by(j);
        k13=k13-bx(j)*bz(j);
        k21=k21-by(j)*bx(j);
        k22=k22+b2(j)-by(j)*by(j);
        k23=k23-by(j)*bz(j);
        k31=k31-bz(j)*bx(j);
        k32=k32-bz(j)*by(j);
        k33=k33+b2(j)-bz(j)*bz(j);
        end
        k11=k11/ncase;
        k12=k12/ncase;
        k13=k13/ncase;
        k21=k21/ncase;
        k22=k22/ncase;
        k23=k23/ncase;
        k31=k31/ncase;
        k32=k32/ncase;
        k33=k33/ncase;
 
%  **   determinate of Ko
 dk=k11*k22*k33+k12*k23*k31+k13*k21*k32-(k13*k22*k31+k11*k32*k23+k12*k21*k33);

 
%  **   cofactor of Ko
        ak11=k22*k33-k23*k32;
        ak12=-(k21*k33-k23*k31);
        ak13=k21*k32-k22*k31;
        ak21=-(k12*k33-k13*k32);
        ak22=(k11*k33-k13*k31);
        ak23=-(k11*k32-k12*k31);
        ak31=k12*k23-k13*k22;
        ak32=-(k11*k23-k13*k21);
        ak33=k11*k22-k12*k21;
 
%  **   calculate inverse Ko
        ik11=ak11/dk;
        ik12=ak21/dk;
        ik13=ak31/dk;
        ik21=ak12/dk;
        ik22=ak22/dk;
        ik23=ak32/dk;
        ik31=ak13/dk;
        ik32=ak23/dk;
        ik33=ak33/dk;
 
%  **   deHoffmann Teller velocity
 
        vxh=ik11*kvx+ik12*kvy+ik13*kvz;
        vyh=ik21*kvx+ik22*kvy+ik23*kvz;
        vzh=ik31*kvx+ik32*kvy+ik33*kvz;
 
 		vht = [vxh vyh vzh]; %km/s
          

         