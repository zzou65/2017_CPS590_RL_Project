%  material derivative of 1-forms

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all
NREFS =4;  % Number of unifrom red refinements
JIG =1;
nsteps=5;

c_Dt = 1;
c_ID = 1;
c_DI = c_ID;

T_Handle=@(t,varargin) cos(2*pi*t);%1/(1+t);%
DT_Handle=@(t,varargin) -2*pi*sin(2*pi*t);%-1/(1+t)^2;%

H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
D11H1=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
D11H2=@(x)(1-x(:,2).^2).*(-2);
D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
D22H2=@(x)(-2).*(1-x(:,1).^2);

H1=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
D1H1=@(x)(-2*x(:,2)).*(-2*x(:,1));
D2H1=@(x)-2*(1-x(:,1).^2);
D11H1=@(x)(-2*x(:,2)).*(-2);
D12H1=@(x)(-2)*(-2*x(:,1));
D22H1=@(x)zeros(size(x,1),1);
H2=@(x) (1-x(:,2).^2).*(2*x(:,1));
D1H2=@(x)(1-x(:,2).^2)*2;
D2H2=@(x)(-2*x(:,2)).*(2*x(:,1));
D11H2=@(x)zeros(size(x,1),1);
D12H2=@(x)(-2*x(:,2))*2;
D22H2=@(x)(-2).*(2*x(:,1));

v1=0.33;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[ H1(x) H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];

%  -v x curl h
IDH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
    DIR1(x).*(D1H2(x)-D2H1(x))];

% grad v h
DIH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
    D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

F_HANDLE=@(x,flag,t,varargin)...
    DtH_HANDLE(x,flag,t) + ...
    IDH_HANDLE(x,flag,t) + ...
    DIH_HANDLE(x,flag,t);

t1=0.3;

% Load mesh from file

% initial mesh
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');

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

directions = Dir_Handle(NewMesh.Coordinates);

M = assemMat_W1F(Mesh,@MASS_W1F,@(x,varargin)1, P3O3());
 
% L = assemLoad_W1F(NewMesh,P7O6(),DtH_HANDLE,t1);
% Dt_h = M\L;
% L = assemLoad_W1F(NewMesh,P7O6(),DIH_HANDLE,t1);
% DI_h=M\L;
% L = assemLoad_W1F(NewMesh,P7O6(),IDH_HANDLE,t1);
% ID_h=M\L;
% L = assemLoad_W1F(NewMesh,P7O6(),H_Handle,t1);
% h_new = M\L;

Dt_h = assemCochain_1f(NewMesh,DtH_HANDLE,gauleg(0,1,10),t1);
DI_h=assemCochain_1f(NewMesh,DIH_HANDLE,gauleg(0,1,10),t1);
ID_h=assemCochain_1f(NewMesh,IDH_HANDLE,gauleg(0,1,10),t1);
h_new = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),t1);

% topological derivatives
TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

% contraction operators
ContrOne=assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms
V=Dir_Handle(NewMesh.Coordinates);
ContrTwo=assemMat_Contr2f(NewMesh,V);               % contraction of two forms

ID = ContrTwo*TopRot;              % -v x curl u geom.
DI = TopGrad*ContrOne;             % grad(v.u) geom.
Err=zeros(nsteps,7);

for k=1:nsteps;

    h = get_MeshWidth(NewMesh);
    dt = 2^-k;
    t0 = t1-dt;

    %     Lh_old = assemLoad_W1F(NewMesh,P7O6(),H_Handle,t0);
    %     h_old = M\Lh_old;

    h_old =  assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),t0);

    pulled_back_vertices=trace_vertices(NewMesh,-dt*directions);
    A = assemMat_SemiLag_W1F(NewMesh, pulled_back_vertices);
    
    pulled_back_bcenters=trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-dt);
    B = assemMat_SemiLagQuad_W1F_bary(NewMesh,  pulled_back_bcenters);    
    B = M\B;
    
    MDh = (h_new-(A*h_old))/dt;
    MDhq = (h_new-(B*h_old))/dt;
    Dth = (h_new-h_old)/dt;
    Lie = (ID+DI)*h_new;

    MDh2 = (h_new-h_old)/dt+Lie;

    Err(k,1) = L2Err_W1F(Mesh,MDh,P7O6(),F_HANDLE,0,t1);
    Err(k,2) = L2Err_W1F(Mesh,MDhq,P7O6(),F_HANDLE,0,t1);
    Err(k,3) = L2Err_W1F(Mesh,MDh2,P7O6(),F_HANDLE,0,t1);
    Err(k,4) = L2Err_W1F(Mesh,Dt_h+DI_h+ID_h,P7O6(),F_HANDLE,0,t1);
    Err(k,5) = L2Err_W1F(Mesh,ID*h_new,P7O6(),IDH_HANDLE,0,t1);
    Err(k,6) = L2Err_W1F(Mesh,DI*h_new,P7O6(),DIH_HANDLE,0,t1);
    Err(k,7) = L2Err_W1F(Mesh,Dth,P7O6(),DtH_HANDLE,0,t1)
end
h = get_MeshWidth(NewMesh);
fig = figure('Name','Error material derivative');
plot(2.^-(1:nsteps),Err(:,1),'r.-',...
    2.^-(1:nsteps),Err(:,2),'go-',...
    2.^-(1:nsteps),Err(:,3),'bx-',...
    2.^-(1:nsteps),Err(:,4),'k+-',...
    2.^-(1:nsteps),Err(:,5),'y*-',...
    2.^-(1:nsteps),Err(:,6),'cs-',...
    2.^-(1:nsteps),Err(:,7),'ms-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf \Delta t, h= ',num2str(h),'}']);
ylabel('{\bf Error}');
legend('Location','Northwest','Lagrage','LagrageQ','Euler+UpFEM','Interpolation','ID','DI','Euler');
