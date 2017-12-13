%  driver for Semi-Lagrange, explizit, operor split

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS =3;
JIG =1;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL =1
T0 = 0;
T1= 0.5;
diff =1;

v1=sqrt(2);
v2=sqrt(2);
Dir_Handle=@(x,t)straight_flow(x);%circ_flow(x);%
F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),2);
H_Handle=@(x,flag,t,varargin)...
    [pulse_2D(x) pulse_2D(x)];
% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

% refine
for j=1:NREFS
     Mesh=refine_REG(Mesh);
     % Mesh preprocessing
     
    switch (JIG)
        case 1
            NewMesh = Mesh;
        case 2
            Loc = get_BdEdges(Mesh);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            FixedPos = zeros(size(Mesh.Coordinates,1),1);
            FixedPos(Loc) = 1;
            NewMesh = jiggle(Mesh,FixedPos);
        case 3
            Loc = get_BdEdges(Mesh);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            FixedPos = zeros(size(Mesh.Coordinates,1),1);
            FixedPos(Loc) = 1;
            NewMesh = smooth(Mesh,FixedPos);
    end
    
end
mw = get_MeshWidth(NewMesh);

h = CFL*mw/norm([v1,v2]);
nsteps = ceil((T1-T0)/h);
h = h-(T0+nsteps*h-T1)/nsteps;
steps = nsteps;

time=zeros(nsteps+1,1);
div=zeros(nsteps+1,3);

NewMesh=add_Edge2Elem(NewMesh);
NewMesh=add_Patches(NewMesh);

% plot_Mesh(NewMesh,'petas');
nEdges= size(NewMesh.Edges,1);

% Assemble Curl-curl matrix, MASS matrix and load vector
C = assemMat_W1F(NewMesh,@STIMA_Curl_W1F,SIGMA_HANDLE,P7O6());
M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P7O6());

% topological derivatives
TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
TopRot=assemMat_TopRot(NewMesh);     % topological Rotation

% contraction operators
ContrOne=assemMat_Contr1f(NewMesh,@(x)Dir_Handle(x));  % contraction of one forms
V=Dir_Handle(NewMesh.Coordinates);
%ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms
ContrTwo=assemMat_Contr2f_fast(NewMesh,Dir_Handle);   % contraction of two forms

ID = M*ContrTwo*TopRot;              % -v x curl u geom.
DI = M*TopGrad*ContrOne;             % grad(v.u) geom.

% pullback of edges
directions = Dir_Handle(NewMesh.Coordinates,0);
pbV = trace_vertices(NewMesh,-h*directions);

pMesh.Coordinates=pbV(:,[1,2]);
pMesh.Elements=NewMesh.Elements;
pMesh = add_Edges(pMesh);
Loc = get_BdEdges(NewMesh);
pMesh.BdFlags = zeros(size(pMesh.Edges,1),1);
pMesh.BdFlags(Loc) = -1;
Meshcell=cell(2);
Meshcell{1}=NewMesh;
Meshcell{2}=pMesh;
plot_Meshes(Meshcell,'as');
%filename=['mesh',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw*steps)),'.eps'];
%print('-depsc', filename);
P = assemMat_SemiLag_W1F(NewMesh, pbV);

% stiffness matrices
ASL = M + h*diff*C;
AIE = M + h*diff*C + h * (ID+DI);
BIE = M ;
AOS = M + h*diff*C;
BOS = M - h*(ID+DI);
%BOS = M + h*(ID');

% timestepping
H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
H_init = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
HSL_old = H_init;
HIE_old = HSL_old;
HOS_old = HSL_old;

mSL =(HSL_old'*M*HSL_old);
mIE = (HIE_old'*M*HIE_old);
mOS = (HOS_old'*M*HOS_old);
div(1,1) =1;
div(1,2) = 1;
div(1,3) = 1;
time(1)=0;

for i = 1:nsteps
    [i i*h nsteps]
    time(i+1)=time(i)+h;

    [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h);
    HIE_new = HSL_new;
    HOS_new = HSL_new;

    %L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h);
    L = assemCochain_1f(NewMesh, F_HANDLE, gauleg(0,1,10),(i)*h);

    LSL = h*M*L+M*P*HSL_old-ASL*HSL_new;
    LIE = h*M*L+BIE*HIE_old-AIE*HIE_new;
    LOS = h*M*L+BOS*HOS_old-AOS*HOS_new;

    HSL_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSL(FreeDofs);
    HIE_new(FreeDofs) = AIE(FreeDofs,FreeDofs)\LIE(FreeDofs);
    HOS_new(FreeDofs) = AOS(FreeDofs,FreeDofs)\LOS(FreeDofs);

    % Update vectors
    HSL_old=HSL_new;
    HIE_old=HIE_new;
    HOS_old=HOS_new;
 
    div(i+1,1) = (HSL_old'*M*HSL_old)/mSL;
    div(i+1,2) = (HIE_old'*M*HIE_old)/mIE;
    div(i+1,3) = (HOS_old'*M*HOS_old)/mOS;
    
end % timestep


figure;
plot(time',div(:,1),'bx-', time',div(:,2),'ro-',time',div(:,3),'gd-');
grid('on');
set(gca,'YScale','log');
xlabel(['\bf t, h=',sprintf('%1.3f',mw),...
 ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps),...
 ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw*steps))]);
ylabel('\bf discrete mass');
legend('location','Northeast','Semi-Lagrangian','Implicit Euler','Operator Splitting');
%filename=['weakdiv',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw*steps)),'.eps'];
%print('-depsc', filename);

%points=find(round(Dummy*100)~=0)