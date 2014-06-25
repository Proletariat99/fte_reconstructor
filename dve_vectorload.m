disp(['Phi at the start of CB is ', num2str(handles.phi,'%3.1f')])
disp(['Theta at the start of CB is ', num2str(handles.theta,'%3.1f')])
gbg_BX=handles.a.data(:,2);
gbg_BY=handles.a.data(:,3);
gbg_BZ=handles.a.data(:,4);
handles.gbg_bc=[gbg_BX(handles.i1:handles.i2) gbg_BY(handles.i1:handles.i2) gbg_BZ(handles.i1:handles.i2)];
handles.theta10 = handles.theta + 1*handles.etxt_ff;
handles.theta20 = handles.theta + 2*handles.etxt_ff;
handles.theta30 = handles.theta + 3*handles.etxt_ff;
handles.phi10 = handles.phi + 1*handles.etxt_ff;
handles.phi20 = handles.phi + 2*handles.etxt_ff;
handles.phi30 = handles.phi + 3*handles.etxt_ff;
%--- xs (formerly zs) components by angle
%--- +10 theta
xx10t = sind(handles.theta10)*cosd(handles.phi);
xy10t = sind(handles.theta10)*sind(handles.phi);
xz10t = cosd(handles.theta10);
%--- +20 theta
xx20t = sind(handles.theta20)*cosd(handles.phi);
xy20t = sind(handles.theta20)*sind(handles.phi);
xz20t = cosd(handles.theta20);
%--- +30 theta
xx30t = sind(handles.theta30)*cosd(handles.phi);
xy30t = sind(handles.theta30)*sind(handles.phi);
xz30t = cosd(handles.theta30);
%--- +10 phi
xx10p = sind(handles.theta)*cosd(handles.phi20);
xy10p = sind(handles.theta)*sind(handles.phi20);
xz10p = cosd(handles.theta);
%--- +20 phi
xx20p = sind(handles.theta)*cosd(handles.phi20);
xy20p = sind(handles.theta)*sind(handles.phi20);
xz20p = cosd(handles.theta);
%--- +30 phi
xx30p = sind(handles.theta)*cosd(handles.phi30);
xy30p = sind(handles.theta)*sind(handles.phi30);
xz30p = cosd(handles.theta);
%--- ys (formerly xs) components by angle (theta)
yx10t = cosd(handles.theta10)*cosd(handles.phi);
yy10t = cosd(handles.theta10)*sind(handles.phi);
yz10t = -sind(handles.theta10);
yx20t = cosd(handles.theta20)*cosd(handles.phi);
yy20t = cosd(handles.theta20)*sind(handles.phi);
yz20t = -sind(handles.theta20);
yx30t = cosd(handles.theta30)*cosd(handles.phi);
yy30t = cosd(handles.theta30)*sind(handles.phi);
yz30t = -sind(handles.theta30);
%--- ys (formerly xs) components by angle (phi)
yx10p = cosd(handles.theta)*cosd(handles.phi10);
yy10p = cosd(handles.theta)*sind(handles.phi10);
yz10p = -sind(handles.theta);
yx20p = cosd(handles.theta)*cosd(handles.phi20);
yy20p = cosd(handles.theta)*sind(handles.phi20);
yz20p = -sind(handles.theta);
yx30p = cosd(handles.theta)*cosd(handles.phi30);
yy30p = cosd(handles.theta)*sind(handles.phi30);
yz30p = -sind(handles.theta);
%--- zs by angle (do we need this?)
% +10 theta xs (gbg = garbage)
gbg_xs10t = [xx10t xy10t xz10t];
gbg_ys10t = [yx10t yy10t yz10t];
gbg_zs10t = cross(gbg_xs10t,gbg_ys10t);
% +20 theta zs (garbage)
gbg_xs20t = [xx20t xy20t xz20t];
gbg_ys20t = [yx20t yy20t yz20t];
gbg_zs20t = cross(gbg_xs20t,gbg_ys20t);
% + 30 theta zs (garbage)
gbg_xs30t = [xx30t xy30t xz30t];
gbg_ys30t = [yx30t yy30t yz30t];
gbg_zs30t = cross(gbg_xs30t,gbg_ys30t);
% +10 phi zs (garbage)
gbg_xs10p = [xx10p xy10p xz10p];
gbg_ys10p = [yx10p yy10p yz10p];
gbg_zs10p = cross(gbg_xs10p,gbg_ys10p);
% +20 phi zs (garbage)
gbg_xs20p = [xx20p xy20p xz20p];
gbg_ys20p = [yx20p yy20p yz20p];
gbg_zs20p = cross(gbg_xs20p,gbg_ys20p);
% + 30 phi zs (garbage)
gbg_xs30p = [xx30p xy30p xz30p];
gbg_ys30p = [yx30p yy30p yz30p];
gbg_zs30p = cross(gbg_xs30p,gbg_ys30p);
% x vectors
handles.bxt1 = (handles.bc.*handles.nb)*gbg_xs10t';
handles.bxt2 = (handles.bc.*handles.nb)*gbg_xs20t';
handles.bxt3 = (handles.bc.*handles.nb)*gbg_xs30t';
handles.bxp1 = handles.bc.*handles.nb*gbg_xs10p';
handles.bxp2 = handles.bc.*handles.nb*gbg_xs20p';
handles.bxp3 = handles.bc.*handles.nb*gbg_xs30p';
% y vectors
handles.byt1 = handles.bc.*handles.nb*gbg_ys10t';
handles.byt2 = handles.bc.*handles.nb*gbg_ys20t';
handles.byt3 = handles.bc.*handles.nb*gbg_ys30t';
handles.byp1 = handles.bc.*handles.nb*gbg_ys10p';
handles.byp2 = handles.bc.*handles.nb*gbg_ys20p';
handles.byp3 = handles.bc.*handles.nb*gbg_ys30p';

%-----------------NEGATIVE ANGLE ADJUSTMENTS---------%
guidata(hObject,handles);                            %
% negative angles                                    %
%--- theta                                           %
ntheta10 = handles.theta - 1*handles.etxt_ff;        %
ntheta20 = handles.theta - 2*handles.etxt_ff;        %
ntheta30 = handles.theta - 3*handles.etxt_ff;        %
%--- phi                                             %
nphi10 = handles.phi - 1*handles.etxt_ff;            %
nphi20 = handles.phi - 2*handles.etxt_ff;            %
nphi30 = handles.phi - 3*handles.etxt_ff;            %
handles.ntheta10=ntheta10;                           %
handles.ntheta20=ntheta20;                           %
handles.ntheta30=ntheta30;                           %
handles.nphi10=nphi10;                               %
handles.nphi20=nphi20;                               %
handles.nphi30=nphi30;      % lazy man's find/replace%
guidata(hObject,handles);       % update handles.
%--- negative xs by angle                            %
%--- xs (formerly zs) (theta)                        %
nxx10t = sind(ntheta10)*cosd(handles.phi);           %
nxy10t = sind(ntheta10)*sind(handles.phi);           %
nxz10t = cosd(ntheta10);                             %
nxx20t = sind(ntheta20)*cosd(handles.phi);           %
nxy20t = sind(ntheta20)*sind(handles.phi);           %
nxz20t = cosd(ntheta20);                             %
nxx30t = sind(ntheta30)*cosd(handles.phi);           %
nxy30t = sind(ntheta30)*sind(handles.phi);           %
nxz30t = cosd(ntheta30);                             %
%--- xs (formerly zs) (phi)
nxx10p = sind(handles.theta)*cosd(nphi10);           %
nxy10p = sind(handles.theta)*sind(nphi10);           %
nxz10p = cosd(handles.theta);                        %
nxx20p = sind(handles.theta)*cosd(nphi20);           %
nxy20p = sind(handles.theta)*sind(nphi20);           %
nxz20p = cosd(handles.theta);                        %
nxx30p = sind(handles.theta)*cosd(nphi30);           %
nxy30p = sind(handles.theta)*sind(nphi30);           %
nxz30p = cosd(handles.theta);                        %
%--- ys (formerly xs) components by angle (theta)
nyx10t = cosd(handles.ntheta10)*cosd(handles.phi);
nyy10t = cosd(handles.ntheta10)*sind(handles.phi);
nyz10t = -sind(handles.ntheta10);
nyx20t = cosd(handles.ntheta20)*cosd(handles.phi);
nyy20t = cosd(handles.ntheta20)*sind(handles.phi);
nyz20t = -sind(handles.ntheta20);
nyx30t = cosd(handles.ntheta30)*cosd(handles.phi);
nyy30t = cosd(handles.ntheta30)*sind(handles.phi);
nyz30t = -sind(handles.ntheta30);
%--- ys (formerly xs) component by angle (nphi)
nyx10p = cosd(handles.theta)*cosd(handles.nphi10);
nyy10p = cosd(handles.theta)*sind(handles.nphi10);
nyz10p = -sind(handles.theta);
nyx20p = cosd(handles.theta)*cosd(handles.nphi20);
nyy20p = cosd(handles.theta)*sind(handles.nphi20);
nyz20p = -sind(handles.theta);
nyx30p = cosd(handles.theta)*cosd(handles.nphi30);
nyy30p = cosd(handles.theta)*sind(handles.nphi30);
nyz30p = -sind(handles.theta);
guidata(hObject,handles);                                        %
% -10 theta zs (garbage)                             %
gbg_nys10t = [nyx10t nyy10t nyz10t];                 %
gbg_nxs10t = [nxx10t nxy10t nxz10t];                 %
gbg_nzs10t = cross(gbg_nxs10t,gbg_nys10t);           %
% -20 theta zs (garbage)                             %
gbg_nys20t = [nyx20t nyy20t nyz20t];                 %
gbg_nxs20t = [nxx20t nxy20t nxz20t];                 %
gbg_nzs20t = cross(gbg_nxs20t,gbg_nys20t);            %
% -30 theta zs (garbage)                             %
gbg_nys30t = [nyx30t nyy30t nyz30t];                 %
gbg_nxs30t = [nxx30t nxy30t nxz30t];                 %
gbg_nzs30t = cross(gbg_nxs30t,gbg_nys30t);            %
% -10 phi zs (garbage)                               %
gbg_nys10p = [nyx10p nyy10p nyz10p];                 %
gbg_nxs10p = [nxx10p nxy10p nxz10p];                 %
gbg_nzs10p = cross(gbg_nxs10p,gbg_nys10p);           %
% -20 phi zs (garbage)                               %
gbg_nys20p = [nyx20p nyy20p nyz20p];                 %
gbg_nxs20p = [nxx20p nxy20p nxz20p];                 %
gbg_nzs20p = cross(gbg_nxs20p,gbg_nys20p);          %
% -30 phi zs (garbage)                               %
gbg_nys30p = [nyx30p nyy30p nyz30p];                 %
gbg_nxs30p = [nxx30p nxy30p nxz30p];                 %
gbg_nzs30p = cross(gbg_nxs30p,gbg_nys30p);           %
% negative x (nx) vectors                            %
handles.nbxt1 = handles.bc.*handles.nb*gbg_nxs10t';  %
handles.nbxt2 = handles.bc.*handles.nb*gbg_nxs20t';  %
handles.nbxt3 = handles.bc.*handles.nb*gbg_nxs30t';  %
handles.nbxp1 = handles.bc.*handles.nb*gbg_nxs10p';  %
handles.nbxp2 = handles.bc.*handles.nb*gbg_nxs20p';  %
handles.nbxp3 = handles.bc.*handles.nb*gbg_nxs30p';  %
% negative y (ny) vectors                            %
handles.nbyt1 = handles.bc.*handles.nb*gbg_nys10t';  %
handles.nbyt2 = handles.bc.*handles.nb*gbg_nys20t';  %
handles.nbyt3 = handles.bc.*handles.nb*gbg_nys30t';  %
handles.nbyp1 = handles.bc.*handles.nb*gbg_nys10p';  %
handles.nbyp2 = handles.bc.*handles.nb*gbg_nys20p';  %
handles.nbyp3 = handles.bc.*handles.nb*gbg_nys30p';  %
%-------END NEGATIVE VECTORS-------------------------%
guidata(hObject,handles);
