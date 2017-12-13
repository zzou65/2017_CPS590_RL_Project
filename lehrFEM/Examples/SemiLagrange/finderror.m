%  driver for Semi-Lagrange, explizit, operor split

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS =1;
JIG =1;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL =10^-6
T0 = 0;
T1= 2;

H1 = @(x) 4 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
D1H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
D2H1=@(x) 4 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
D11H1=@(x) -16 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
D12H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
D22H1=@(x) -24 * (1-x(:,1).^2).^2 .*x(:,2);
H2=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
D1H2=@(x) -4 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
D2H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
D11H2=@(x) 24 * x(:,1).*(1-x(:,2).^2).^2;
D12H2=@(x) 16 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
D22H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
T_Handle=@(t,varargin) cos(2*pi*t);%1/(1+t);%

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[ H1(x) H2(x)];
F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),2);

% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

Mesh=refine_REG(Mesh);     

NewMesh = load_Mesh('Coord.dat','Elem.dat');% initial mesh
NewMesh = add_Edges(NewMesh);
Loc = get_BdEdges(NewMesh);
NewMesh.BdFlags = zeros(size(NewMesh.Edges,1),1);
NewMesh.BdFlags(Loc) = -1;

%NewMesh.Coordinates(1:21,:)=Mesh.Coordinates(1:21,:)
%NewMesh.Coordinates(23:end,:)=Mesh.Coordinates(23:end,:)
NewMesh.Coordinates=Mesh.Coordinates

range=[0:0.0001:0.025];
range=[0:0.001:0.05];
pbV=cell(size(range,2),1);
P=cell(size(range,2),1);
M=cell(size(range,2),1);

div=zeros(size(range,2),1);

for r=1:size(range,2)
    [r,size(range,2)]
    NewMesh.Coordinates(15,:)=Mesh.Coordinates(15,:)-[10^-12,-range(r)];%round(100*NewMesh.Coordinates(22,:))/100;

    mw = get_MeshWidth(NewMesh);

    h = CFL*mw/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps;
    steps = nsteps;

    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % plot_Mesh(NewMesh,'petas');
    nEdges= size(NewMesh.Edges,1);

    % Assemble Curl-curl matrix, MASS matrix and load vector
    M{r} = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P7O6());

    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient

    % pullback of edges
    directions = Dir_Handle(NewMesh.Coordinates,0);
    pbV{r} = pullback_Vertices(NewMesh,h*directions);

    pMesh.Coordinates=pbV{r}(:,[1,2]);
    pMesh.Elements=NewMesh.Elements;
    pMesh = add_Edges(pMesh);
    Loc = get_BdEdges(NewMesh);
    pMesh.BdFlags = zeros(size(pMesh.Edges,1),1);
    pMesh.BdFlags(Loc) = -1;
    plot_pbMesh(NewMesh,pMesh,pMesh,'pa');

    P{r} = pullback_Edges(NewMesh, pbV{r});

    %test_pbE(NewMesh,pMesh,P);

    % timestepping
    H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
    HSL_old = M{r}\H_init;

    [Dummy,FD_LFE] = assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
    Dummy(FD_LFE) =TopGrad(:,FD_LFE)'*M{r}*HSL_old;
    %Dummy(FD_LFE) =TopGrad(:,FD_LFE)'*P*M*HSL_old-TopGrad(:,FD_LFE)'*M*HSL_old;
   % plot_LFE(Dummy,NewMesh);colorbar;
    %div(r) = norm(TopGrad(:,FD_LFE)'*M{r}*HSL_old);
    time(1)=0;

    time(2)=time(1)+h;

    [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h);

    L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h);

    LSL = h*L-P{r}'*M{r}*HSL_old-M{r}*HSL_new;

    HSL_new(FreeDofs) = M{r}(FreeDofs,FreeDofs)\LSL(FreeDofs);

    HSL_old=HSL_new;
    Dummy(FD_LFE)=TopGrad(:,FD_LFE)'*M{r}*HSL_old;
    %plot_LFE(Dummy,NewMesh); colorbar;
    %plot_Norm_W1F(M{r}*HSL_old,NewMesh); colorbar;
    div(r) = norm(TopGrad(:,FD_LFE)'*M{r}*HSL_old);

end % timestep
figure;
plot(range,div);