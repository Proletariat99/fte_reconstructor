function [ehts,ecs,R,slope,vht]=HTcoef(bc,vp,ndata)
e0=1e-03;

vht=dht(bc,vp,ndata);
%  **   Eht=-Vht x B ,Ec=-U x B
for j=1:ndata
    ecx(j)=e0*(bc(j,2)*vp(j,3)-bc(j,3)*vp(j,2));
    ecy(j)=e0*(bc(j,3)*vp(j,1)-bc(j,1)*vp(j,3));
    ecz(j)=e0*(bc(j,1)*vp(j,2)-bc(j,2)*vp(j,1));
    ehtx(j)=e0*(bc(j,2)*vht(1,3)-bc(j,3)*vht(1,2));
    ehty(j)=e0*(bc(j,3)*vht(1,1)-bc(j,1)*vht(1,3));
    ehtz(j)=e0*(bc(j,1)*vht(1,2)-bc(j,2)*vht(1,1));
end

ecs=[ecx' ecy' ecz'];
ehts=[ehtx' ehty' ehtz'];

%  **   Correlation Coefficient
Eht=zeros(1,ndata*3);
for  j=1:ndata
    Eht(j)        =Eht(j)        +ehtx(j);
    Eht(j+ndata)  =Eht(j+ndata)  +ehty(j);
    Eht(j+ndata*2)=Eht(j+ndata*2)+ehtz(j);
end
Ec=zeros(1,ndata*3);
for  j=1:ndata
    Ec(j)        =Ec(j)        +ecx(j);
    Ec(j+ndata)  =Ec(j+ndata)  +ecy(j);
    Ec(j+ndata*2)=Ec(j+ndata*2)+ecz(j);
end
ttemp=corrcoef(Eht,Ec);
R=ttemp(1,2);
[p,s]=polyfit(Eht,Ec,1);
yy=polyval(p,Eht);
slope=p(1);
