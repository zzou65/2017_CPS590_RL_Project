%  Test pullback of vertieces

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

close all;
clear Mesh;
NREFS = 2;  % Number of unifrom red refinements
JIG =1;
%Dir_Handle=@(x)-[0.1*ones(size(x,1),1), 0.25*ones(size(x,1),1)];
%Dir_Handle=@(x)-[0.1*(1+x(:,1)+x(:,2)+x(:,1).*x(:,2)).*ones(size(x,1),1), zeros*(1+(x(:,1)+x(:,2)+x(:,1).*x(:,2))).*ones(size(x,1),1)];
%Dir_Handle=@(x,varargin)10^-2*straight_flow(x);
%Dir_Handle=@(x)-[0.1*(1+x(:,1)+x(:,2)+x(:,1).*x(:,2)).*ones(size(x,1),1), 0.2*(1+(x(:,1)+x(:,2)+x(:,1).*x(:,2))).*ones(size(x,1),1)];
Dir_Handle=@(x,varargin) -(1-x(:,2).^2).*(1-x(:,1).^2)*[0.4 -0.2];

G_Handle=@(x,flag,t,varargin) cos(2*pi*t)*(1+x(:,1)+x(:,2));
DtG_Handle=@(x,flag,t,varargin) -2*pi*sin(2*pi*t)*(1+x(:,1)+x(:,2));
LtG_Handle=@(x,flag,t,varargin) cos(2*pi*t)*0.2*(-1+x(:,2).^2).*(1-x(:,1).^2).*(1+x(:,1)+x(:,2));

% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

% refine
for i=1:NREFS
    Mesh=refine_REG(Mesh);
end

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


NewMesh=add_Edge2Elem(NewMesh);
NewMesh=add_Patches(NewMesh);

plot_Mesh2(NewMesh,'petas');

%plot velocities
figure;
M = assemMat_W1F(Mesh,@MASS_W1F,@(x,varargin)1, P3O3());
Dir = assemLoad_W1F(Mesh,P7O6(),Dir_Handle,0);
plot_W1F(M\Dir,NewMesh);

CFL =[0.4];
T0 = 0;
T1= 2;
v1=1/sqrt(2);
v2=1/sqrt(2);
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(881-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];

mw = get_MeshWidth(NewMesh);
h = CFL*mw/norm([v1,v2]);
nsteps = ceil((T1-T0)/h);
h = h-(T0+nsteps*h-T1)/nsteps;

% pullback of vertices
directions=Dir_Handle(NewMesh.Coordinates);
pulled_back_vertices=trace_vertices(NewMesh,h*directions);
pulled_back_bcenters=trace_bcenters(NewMesh,0*directions);
B1q = assemMat_SemiLagQuad_W1F_bary(NewMesh, pulled_back_bcenters,JDir_Handle,0);

M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P1O2());
%A =assemMat_pbV(Mesh,pulled_back_vertices);

% [U,D]=eig(full(A));
% for i =1:size(Mesh.Coordinates,1)
%     plot_LFE(U(:,i),NewMesh); colorbar
% end
plot_transportMesh(NewMesh,pulled_back_vertices,'tas');
% plot_pbMesh(NewMesh,pMesh{j},'a');
%plot_Mesh2(pMesh,'petas');

%pbE=assemMat_pbE(NewMesh, pulled_back_vertices);
%pbE=pullback_Edges(NewMesh, pulled_back_vertices);
