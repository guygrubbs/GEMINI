function [z,x,y]=geog2UEN(alt,glon,glat,thetactr,phictr)

%DISTANCE EW AND NS (FROM ENU (or UEN in our case - cyclic permuted) COORD. SYSTEM) NEED TO BE CONVERTED TO DIPOLE SPHERICAL AND THEN
%GLAT/GLONG - BASICALLY HERE WE ARE MAPPING THE CARTESIAN GRID ONTO THE
%SURFACE OF A SPHERE THEN CONVERTING TO GEOGRAPHIC.

%UPWARD DISTANCE
Re=6370e3;
z=alt;

%Convert to geomganetic coordinates
[theta,phi]=geog2geomag(glat,glon);


%Convert ot northward distance in meters
gamma2=theta-thetactr;    %southward magnetic angular distance
gamma2=-1*gamma2;    %convert to northward angular distance
y=gamma2*Re;

%{
%Northward angular distance
gamma2=y/Re;    %must retain the sign of x3
theta=thetactr-gamma2;   %minus because distance north is against theta's direction
theta=reshape(theta,[1,1,lx3]);
theta=repmat(theta,[lx1,lx2,1]);
%}


gamma1=phi-phictr;    %eastward angular distance
x=Re*sin(thetactr)*gamma1;

%{
%Eastward angular distance
%gamma1=x/Re;     %must retain the sign of x2
gamma1=x/Re/sin(thetactr);     %must retain the sign of x2, just use theta of center of grid
phi=phictr+gamma1;
phi=reshape(phi,[1,lx2,1]);
phi=repmat(phi,[lx1,1,lx3]);
%}
end
