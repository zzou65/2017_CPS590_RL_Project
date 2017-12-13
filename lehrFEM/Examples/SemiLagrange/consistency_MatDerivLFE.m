%  material derivative of 0-forms, consistency of temporal derivative ,
%  material derivative ...  mesh and time refinement

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS = 3;  % Number of unifrom red refinements
JIG =1;
Ndt=16;     % Number of delta 
offset=5;

T_Handle=@(t,varargin) t^3;
DT_Handle=@(t,varargin) 3*t^2;
v1=0.4;
v2=0.2;

Dir_Handle=@(x,varargin) -(1-x(:,2).^2).*(1-x(:,1).^2)*[v1 v2];

G_Handle=@(x,flag,t,varargin) T_Handle(t)*(1+x(:,1).^2+x(:,2).^2);
DtG_Handle=@(x,flag,t,varargin) DT_Handle(t)*(1+x(:,1).^2+x(:,2).^2);
LtG_Handle=@(x,flag,t,varargin) T_Handle(t)*(v1*2*x(:,1)+v2*2*x(:,2)).*(-1+x(:,2).^2).*(1-x(:,1).^2);
F_Handle=@(x,flag,t,varargin) DtG_Handle(x,flag,t)+LtG_Handle(x,flag,t); 

t1=0.24;
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

% topological derivatives
%TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient

% contraction operators
%ContrOne=assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms

M_LFE = assemMat_LFE(NewMesh,@MASS_LFE);
M_0 =  assemMat_Mass0fD(NewMesh);
%  
%standard FEM
ID = assemMat_LFE(NewMesh, @STIMA_Conv_LFE,Dir_Handle,P3O3());
ID = M_LFE\ID;
% upwind FEM
ID_up = assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, Dir_Handle);
%ID_up = ContrOne*TopGrad;

% L1 = assemLoad_LFE(NewMesh,P7O6(),DtG_Handle,t1);
% g1=M_LFE\L1;
% L2 = assemLoad_LFE(NewMesh,P7O6(),LtG_Handle,t1);
% g2=M_LFE\L2;
% Lg_new= assemLoad_LFE(NewMesh,P7O6(),G_Handle,t1);
% g_new=M_LFE\Lg_new;

g1 = assemCochain_0f(NewMesh,DtG_Handle,t1);
g2 = assemCochain_0f(NewMesh,LtG_Handle,t1);
g_new= assemCochain_0f(NewMesh,G_Handle,t1);

Err=zeros(Ndt,6);

for k=(1:Ndt);
    
    dt=2^-(offset+k);
    t0=t1-dt;

    % explizit Euler to trace characteristics
    directions=Dir_Handle(NewMesh.Coordinates);
    pulled_back_vertices=trace_vertices(NewMesh,-dt*directions);
    A =assemMat_SemiLag_LFE(Mesh,pulled_back_vertices);
%     B =assemMat_SemiLagQuad_LFE(Mesh,pulled_back_vertices);
%     A=M_0\B;
    
    %     Lg_old= assemLoad_LFE(NewMesh,P7O6(),G_Handle,t0);
    %     g_old=M_LFE\Lg_old;
    g_old= assemCochain_0f(NewMesh,G_Handle,t0);
   
    MDg = (g_new-A*g_old)/dt;
    Dtg = (g_new-g_old)/dt;
    MDg2 = Dtg+ID*g_new;
    MDg3 = Dtg+ID_up*g_new;
    
    Err(k,1)=L2Err_LFE(NewMesh,MDg,P7O6(),F_Handle,0,t1);
    Err(k,2)=L2Err_LFE(NewMesh,MDg2,P7O6(),F_Handle,0,t1);
    Err(k,3)=L2Err_LFE(NewMesh,MDg3,P7O6(),F_Handle,0,t1);
    Err(k,4)=L2Err_LFE(NewMesh,g1+g2,P7O6(),F_Handle,0,t1);
    Err(k,5)=L2Err_LFE(NewMesh,ID_up*g_new,P7O6(),LtG_Handle,0,t1);
    Err(k,6)=L2Err_LFE(NewMesh,Dtg,P7O6(),DtG_Handle,0,t1)
end
h = get_MeshWidth(NewMesh);
fig = figure('Name','Consistency error material derivative');
plot(2.^-(offset+(1:Ndt)),Err(:,1),'r.-',...
    2.^-(offset+(1:Ndt)),Err(:,2),'bo-',...
    2.^-(offset+(1:Ndt)),Err(:,3),'gx-',...
    2.^-(offset+(1:Ndt)),Err(:,4),'y+-',....
    2.^-(offset+(1:Ndt)),Err(:,5),'c*-',....
    2.^-(offset+(1:Ndt)),Err(:,6),'ks-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t}');
xlabel(['{\bf \Delta t, h= ',num2str(h),'}']);
ylabel('{\bf Error}');
legend('Location','Northwest','Lagrange','Euler+FEM','Euler+UpFEM', 'Interpolation','UpFEM','Euler');

% fig = figure('Name','Error time derivative');
% plot(2.^(offset+(1:Ndt)),Err(:,6),'ro--',...
%     2.^(offset+(1:Ndt)),Err(:,7),'bo--'); grid('on');
% set(gca,'YScale','log','XScale','log');
% xlabel('{\bf \Delta t}');
% ylabel('{\bf Error}');
