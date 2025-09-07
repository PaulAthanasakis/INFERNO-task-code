clear;clc; close all;

%% Define Model Variables
Length=0.3; %Assume total length of canal 30cm
H=0.03; %Width of canal in m

T_Hot = 40; %degC
T_Cold= 10; %degC
T_Env=17; %degC
%Assume water passes through the channel
m= 0.891*0.001; %Pas
k= 0.607; %W/mK
Cp=4180; %J/kgK
den=997; %kg/m^3
a=k/den/Cp;
dP_dx = -2; %pressure difference in channel

%% Define Grid Mesh
dh=1 ;
L=100; 
x=0:dh:10*L; x(end)=10*L; 
y=0:dh:L; y(end)=L; 
nx=length(x);  % number of nodes in x direction
ny=length(y);  % number of nodes in y direction 
dh=Length/(nx-1);
n=nx*ny;       % total number of nodes 

%% Mark Nodes
top=ny+1:ny:n-2*ny+1; bottom=2*ny:ny:n-ny;
left=2:ny-1; right=(n-ny+2):n-1;
corners=[1 ny n-ny+1 n];
inner=1:n; inner([left right top bottom corners])=[];
border=1:n; border(inner)=[];

plate_cover = 0.3; %percentage of plates to total length
plate_length = plate_cover*Length;
HotPlate = [corners(1) top(1:plate_length/dh)]; %nodes corresponding to Hot Plate
ColdPlate = [corners(2) bottom(1:plate_length/dh)]; %nodes corresponding to Cold Plate

%% Peclet Number

v_mean=(-1/12)*(1/m)*(dP_dx)*H^2; 
Re=den*v_mean*10*L/m;
Pr=m*Cp/k;
Pe=Re*Pr;

%% Calculations 

h=H/(ny-1).*y;
B=-(dh)^2*(m/k)*(1/(2*m))*(dP_dx).*(h.^2-2*H.*h+H^2);
B=B'; B=flip(B); 
B=repmat(B,1,nx);
B(HotPlate)=T_Hot; B(ColdPlate)=T_Cold; 
B(left)=T_Env;
B=reshape(B,n,1);

vx=-(1/8)*(dP_dx)*H^2.*(1-(2.*(h/H)-1).^2);
c=dh/a.*vx;
c1=1+c;  c2=4+c; c3=2+c;

% D_minus_ny=[repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) ones(1,n-ColdPlate(end)-ny) 2*ones(1,ny) zeros(1,ny)];
% D_minus_1=[[0 zeros(1,ny-2) 0] repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([0 ones(1,ny-2) 2],1,nx-length(ColdPlate)) 0];
% D_minus_1(1)=[];
% D_0=[[1 ones(1,ny-2) 1] repmat([1 -4*ones(1,ny-2) 1],1,length(ColdPlate)-1) -4*ones(1,n-ColdPlate(end))];
% D_plus_1=[0 [0 zeros(1,ny-2) 0] repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([2 ones(1,ny-2) 0],1,nx-length(ColdPlate))]; D_plus_1(end)=[];
% D_plus_ny=[zeros(1,ny) zeros(1,ny) repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) ones(1,n-ColdPlate(end)-ny)];

D_minus_ny=[ repmat([0 c1(2:end-1) 0],1,length(ColdPlate)-1) repmat(c1,1,nx-length(ColdPlate)-1) c3 zeros(1,ny)];
D_minus_1=[zeros(1,ny-1) repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([0 ones(1,ny-2) 2],1,nx-length(ColdPlate)) 0];
D_0=[ones(1,ny) repmat([1 -c2(2:end-1) 1],1,length(ColdPlate)-1) repmat(-c2,1,nx-length(ColdPlate))];
D_plus_1=[0 zeros(1,ny) repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([2 ones(1,ny-2) 0],1,nx-length(ColdPlate))]; D_plus_1(end)=[];
D_plus_ny=[zeros(1,2*ny) repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) ones(1,n-ColdPlate(end)-ny)];

A=sparse(n,n);
A=spdiags(D_minus_ny',-ny,A);
A=spdiags(D_minus_1',-1,A);
A=spdiags(D_0',0,A);
A=spdiags(D_plus_1',1,A);
A=spdiags(D_plus_ny',ny,A);

T=A\B;
T=reshape(T,ny,nx);
figure(1);contourf(T,"ShowText",true); 
set(gca,'YDir','normal')
title('Temperature Profile within entire channel')
xlabel('Length of channel (cm)')
xticks(floor(nx/10*(0:1:10)))
xticklabels({'0','3','6','9','12','15','18','21','24','27','30'})
ylabel('width of channel (cm)')
yticks(floor(ny/10*(0:1:10)))
yticklabels({'0','0.3','0.6','0.9','1.2','1.5','1.8','2.1','2.4','2.7','3.0'})

% 
% T_graph=T(:,[1 floor(length(ColdPlate)/2) length(ColdPlate)+3 floor(0.8*nx) nx]);
% figure(2); plot(T_graph)
% legend({'entrance','within plates','right after plates','away from plates','end of channel'})
% title('Temperature Profiles at select cross sections')
% xlabel('Cross Section of Sensor Array-nodes in y direction')
% ylabel('Temperature (C)')
% axis([1 ny 10 40])