function  [A BJ EJ] = MaxwellTimeStep(Sigma,A_old,Mesh,t,dt,v,A_Bd_Handle,Load_Handle,VxCurl_Handle,theta)
%  performs timestep and returns A and coupling values BxJ and E J.
%  INPUT:
%   Sigma: n-1 array, n number of gridpoints
%   A_Old: m-1 array, m number of edges
%   t:  total time   
%   dt: timestepsize
%   v:   n-2 array, velocity at gripoints
%   A_Bd_Handle: function handle to handle A on inflow boundary 
%   Load_Handle: function handle to handle righthand side
%   VxCurl_Handle: function handle to handle VxCurl on inflow boundary 
%
%   OutPut:
%   A: m-1 array iterated vectorpotential
%   BJ: n-2 array
%   EJ: n-1 array
%   Example:
%
%   

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

QuadRule=Duffy(TProd(gauleg(0,1,10)));
Zero_Handle=@(x,varargin) zeros(size(x,1), 2);
% sigma-weighted Mass of one forms
M_sigma=assemMat_W1Fdis(Mesh,Sigma);
M_W1F=assemMat_W1Fdis(Mesh,Sigma);

% topological curl operator
TopRot=assemMat_TopRot(Mesh);

% contraction of two forms 
ContrTwo=assemMat_ContrTwo(Mesh,v);

% diagonal mass matrix of two-forms
MassTwo=assemMat_MassTwoD(Mesh);

% Stiffnessmatrix
%S=M_sigma+dt*TopRot'*MassTwo*TopRot+dt*M_sigma*ContrTwo*TopRot;
S=M_sigma+dt*TopRot'*MassTwo*TopRot;

%iterated A
FreeDofs=intersect( find(Mesh.BdFlags ~=-1), find(Mesh.BdFlags ~=-2));
nonInFlow=find(Mesh.BdFlags ~=-1);
% nonOutFlow=find(Mesh.BdFlags ~=-1);

cL=assemCochain_1f(Mesh,Load_Handle,gauleg(0,1,10),t); 
cL=dt*cL;
A=assemCochain_1f(Mesh,A_Bd_Handle,gauleg(0,1,10),t);
A(FreeDofs)=0;
c1=assemCochain_1f(Mesh,VxCurl_Handle,gauleg(0,1,10),t);
c1(nonInFlow)=0;

% % IMLIZIT
% Lf=M_sigma*(dt*cL-dt*c1+A_old)-S*A;
% A(FreeDofs)=S(FreeDofs,FreeDofs)\Lf(FreeDofs);
Lf=M_sigma*(dt*cL-A_old)-S*A;
A(FreeDofs)=S(FreeDofs,FreeDofs)\Lf(FreeDofs);


% % EXLIZIT
% Lf=M_sigma*(dt*cL)+S*A_old;
% A(FreeDofs)=M_sigma(FreeDofs,FreeDofs)\Lf(FreeDofs);

% elctric field
E=(-A+A_old)/dt;

% electric current
J=TopRot'*MassTwo*TopRot*A;

% magnetic field
B=TopRot;

% Energy
EJ=assemWedge_1f1f(Mesh,E,J);

% ?? Momentum??
BJ=assemWedge_2f1f(Mesh,B,J);
return;