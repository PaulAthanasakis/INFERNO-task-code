# INFERNO-task-code
Code for solving a 2d diffusion-convection of heat in a channel 
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
dP_dx = -3.6; %pressure difference in channel

%% Define Grid Mesh
dh=1 ;
L=10; 
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

v_mean=(-1/12)*(dP_dx)*H^2; 
Re=den*v_mean*10*L/m;
Pr=m*Cp/k;
Pe=Re*Pr;

%% Calculations 

B=-(dh)^2*(m/k)*(1/(2*m))*(dP_dx).*(y.^2-2*H.*y+H^2);
B=B'; B=flip(B); 
B=repmat(B,1,nx);
B(HotPlate)=T_Hot; B(ColdPlate)=T_Cold; 
B(left)=T_Env;
B=reshape(B,n,1);

D_minus_ny=[repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) ones(1,n-ColdPlate(end)-ny) 2*ones(1,ny) zeros(1,ny)];
D_minus_1=[[0 zeros(1,ny-2) 0] repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([0 ones(1,ny-2) 2],1,nx-length(ColdPlate)) 0];
D_minus_1(1)=[];
D_0=[[1 ones(1,ny-2) 1] repmat([1 -4*ones(1,ny-2) 1],1,length(ColdPlate)-1) -4*ones(1,n-ColdPlate(end))];
D_plus_1=[0 [0 zeros(1,ny-2) 0] repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) repmat([2 ones(1,ny-2) 0],1,nx-length(ColdPlate))]; D_plus_1(end)=[];
D_plus_ny=[zeros(1,ny) zeros(1,ny) repmat([0 ones(1,ny-2) 0],1,length(ColdPlate)-1) ones(1,n-ColdPlate(end)-ny)];

A=sparse(n,n);
A=spdiags(D_minus_ny',-ny,A);
A=spdiags(D_minus_1',-1,A);
A=spdiags(D_0',0,A);
A=spdiags(D_plus_1',1,A);
A=spdiags(D_plus_ny',ny,A);

T=A\B;
T=reshape(T,ny,nx);
figure(1);contourf(T)
% T_graph=T(:,[1 floor(length(ColdPlate)/2) length(ColdPlate)+3 floor(0.8*nx) nx]);
% figure(2); plot(T_graph')
