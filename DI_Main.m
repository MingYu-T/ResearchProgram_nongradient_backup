%   Displacement Inverter Topology Optimization
%   TO framework    : MMC-Based Method
%   MMC component   : Moving Wide B-Spline with Constrained Ends(fixed width)
%   Algorithm       : GA
%   Note            : Adjust Spline Constrained Ends for different cases!!!
clear;clc;close all
dbstop if error
% - Global Variables
global Ncp N ord E0 nu nele thickness...
    gobalcoords0 edofMat fixeddof loaddof F ndof Kin Kout outdof varfail nelx nely volfrac GenBest_x GenBest_f splineR gen

% - Defining Initial Parameters
%  Ncp      : Number of control points
%  N      : Number of splines
%  splineR        : half spline width  
%  ord            : Order of polynomial spline
%  E0             : Elastic (yield) modulus [Pa]
%  nu             : Poisson's ratio
%  DW             : Design domain width 
%  DH             : Design domain height
%  nelx           : Number of elements on x axis
%  nely           : Number of elements on y axis

Ncp = 5;
N = 6;
splineR = 0.15e-3;
ord = 3;
E0=1e8;
nu=0.3;
nelx = 100;
nely = 50;
volfrac = 0.25;
GenBest_x = [];
GenBest_f = [];
gen = 0;

%  - Initialize Boundary Bondition and Mesh
[nele,thickness,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter;

% - Initial Parameter Limits
% order of var: [spline1(Bx2...Bxn-1 By2...Byn-1 Br1) spline2( ... Br2)... splinen( ... Brn)]
sv_num = (Ncp-2)*2+1; % number of var for 1 spline
v_num = sv_num*N;

lb = zeros(1, v_num);
for i = 1:N
    lb( 1+(i-1)*sv_num:1+(i-1)*sv_num+(Ncp-3) ) = 0; % Bx
    lb( (1+Ncp-2)+(i-1)*sv_num:(1+Ncp-2)+(i-1)*sv_num+(Ncp-3) ) = 0; % By
    lb(sv_num*i) = 0; % Br
end
ub = zeros(1, v_num);
for i = 1:N
    ub( 1+(i-1)*sv_num:1+(i-1)*sv_num+(Ncp-3) ) = 10e-3; % Bx
    ub( (1+Ncp-2)+(i-1)*sv_num:(1+Ncp-2)+(i-1)*sv_num+(Ncp-3) ) = 5e-3; % By
    ub(sv_num*i) = 0.3e-3; % Br
end


% - Topology Synthesis Optimization (Calling for Genetic Algorithm)
tic
rng default
opt_ga = optimoptions(@ga, 'PopulationSize', 150, 'MaxGenerations', 40, ...
   'PlotFcn', {@gaplotdistance, @gaplotbestf}, 'MaxStallGenerations', 250,...
   'OutputFcn',@GetData);
[x, fval, exitflag, output, population, fitness] = ga(@objfunc, v_num, [], [], [], [], lb, ub, [], opt_ga);
toc








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nele,thickness,gobalcoords0,edofMat,fixeddof,loaddof,F,ndof,Kin,Kout,outdof]=initialize2Dinverter
DW = 10e-3; DH = 5e-3;
nelx=100;nely=50;
EW=DW/nelx; EH=DH/nely; thickness=1e-3;
[i0,j0] = meshgrid(0:nelx,0:nely);
gobalcoords0x=i0*EW;
gobalcoords0y=EH*nely-j0*EH;
gobalcoords0xy=[gobalcoords0x(:) gobalcoords0y(:)]';
gobalcoords0=gobalcoords0xy(:);
[il,jl] = meshgrid(0,nely);
loadnid = il*(nely+1)+(nely+1-jl); 
loaddof = 2*loadnid(:)-1; 
force=5;
Kin=500;    % 500 N/m
[ilo,jlo] = meshgrid(nelx,nely);
outnid = ilo*(nely+1)+(nely+1-jlo);  
outdof = 2*outnid(:)-1; 
Kout=100;   % 100 N/m
[iif,jf] = meshgrid(0,0:nely/5);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); 2*fixednid(:)-1];
[iif,jf] = meshgrid(0:nelx,nely);
fixednid = iif*(nely+1)+(nely+1-jf);
fixeddof = [2*fixednid(:); fixeddof];
nele = nelx*nely;
ndof = 2*(nelx+1)*(nely+1);
F = sparse(loaddof,1,force,ndof,1);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);                               
edofVec = 2*nodeids(:)+1;                                                              
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*(nely+1) 2*(nely+1)+1 2*(nely+1)-2 2*(nely+1)-1 -2 -1],nele,1);%%%%b=[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],超神节点提取符写法，阐述其相对单元的1节点x方向自由度的位置
end

function [state,options,optchanged] = GetData(options,state,~)
    global GenBest_x GenBest_f gen
    [best_val,idx] = min(state.Score);
    best_x = state.Population(idx,:);
    GenBest_x = [GenBest_x; best_x];
    GenBest_f = [GenBest_f; best_val];
    TopoCheck(best_x);
    gen = gen+1;
    writematrix(GenBest_x,'b_var.txt','Delimiter','tab')
    optchanged = false;
end