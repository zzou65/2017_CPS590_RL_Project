%  material derivative of 1-forms, consistency,
%  time refinement

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all
NREFS =5;  % Number of unifrom red refinements
JIG =1;

Err=zeros(NREFS,7);
mw=zeros(NREFS,1);
CFL=0.5;

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

% H1=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% D1H1=@(x)(-2*x(:,2)).*(-2*x(:,1));
% D2H1=@(x)-2*(1-x(:,1).^2);
% D11H1=@(x)(-2*x(:,2)).*(-2);
% D12H1=@(x)(-2)*(-2*x(:,1));
% D22H1=@(x)zeros(size(x,1),1);
% H2=@(x) (1-x(:,2).^2).*(2*x(:,1));
% D1H2=@(x)(1-x(:,2).^2)*2;
% D2H2=@(x)(-2*x(:,2)).*(2*x(:,1));
% D11H2=@(x)zeros(size(x,1),1);
% D12H2=@(x)(-2*x(:,2))*2;
% D22H2=@(x)(-2).*(2*x(:,1));

v1=0.66;
v2=0.3;
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

% curl curl h
CCH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22H1(x) + D12H2(x) ...
    D12H1(x)-D11H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];
   
% v * div h
IDH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [DIR1(x).*(D1H1(x)+D2H2(x)) ...
     DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h) 
DIH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
    +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];
 
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
for k=1:NREFS
    Mesh=refine_REG(Mesh);

    
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


    M = assemMat_W1F(Mesh,@MASS_W1F,@(x,varargin)1, P3O3());

    %     L = assemLoad_W1F(NewMesh,P7O6(),DtH_HANDLE,t1);
    %     Dt_h = M\L;
    %     L = assemLoad_W1F(NewMesh,P7O6(),DIH_HANDLE,t1);
    %     DI_h=M\L;
    %     L = assemLoad_W1F(NewMesh,P7O6(),IDH_HANDLE,t1);
    %     ID_h=M\L;
    %     L = assemLoad_W1F(NewMesh,P7O6(),H_Handle,t1);
    %     h_new = M\L;

    Dt_h = assemCochain_1f(NewMesh,DtH_HANDLE,gauleg(0,1,10),t1);
    DI_h=assemCochain_1f(NewMesh,DIH_HANDLE,gauleg(0,1,10),t1);
    ID_h=assemCochain_1f(NewMesh,IDH_HANDLE,gauleg(0,1,10),t1);
    h_new = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),t1);
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

    % contraction operators
    ContrOne=assemMat_Contr1f(NewMesh,@(x)Dir_Handle(x));  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,V);               % contraction of two forms

    ID = ContrTwo*TopRot;              % -v x curl u geom.
    DI = TopGrad*ContrOne;             % grad(v.u) geom.
  
    h = get_MeshWidth(NewMesh);
    dt = CFL*h;
    t0 = t1-dt;

    %Lh_old = assemLoad_W1F(NewMesh,P7O6(),H_Handle,t0);
    %h_old = M\Lh_old;
    h_old =  assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),t0);

    pulled_back_vertices=trace_vertices(NewMesh,dt*V);
    A = assemMat_SemiLag_W1F(NewMesh, pulled_back_vertices);
    
    pulled_back_bcenters=trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,dt);
    B = assemMat_SemiLagQuad_W1F_bary(NewMesh,  pulled_back_bcenters);    
    B = M\B;
   
    MDh = (h_new-(M\(A'*M)*h_old))/dt;
    MDhq = (h_new-((M\B')*h_old))/dt;
    Dth = (h_new-h_old)/dt;
    Lie = M\(ID'*M+DI'*M)*h_new;

    MDh2 = (h_new-h_old)/dt-Lie;

    mw(k)=h;
    Err(k,1) = L2Err_W1F(NewMesh,MDh,P7O6(),F_HANDLE,0,t1);
    Err(k,2) = L2Err_W1F(NewMesh,MDhq,P7O6(),F_HANDLE,0,t1);
    Err(k,3) = L2Err_W1F(NewMesh,MDh2,P7O6(),F_HANDLE,0,t1);
    Err(k,4) = L2Err_W1F(NewMesh,Dt_h+DI_h+ID_h,P7O6(),F_HANDLE,0,t1);
    Err(k,5) = L2Err_W1F(NewMesh,-M\DI'*M*h_new,P7O6(),IDH_HANDLE,0,t1);
    Err(k,6) = L2Err_W1F(NewMesh,-M\ID'*M*h_new,P7O6(),DIH_HANDLE,0,t1);
    Err(k,7) = L2Err_W1F(NewMesh,Dth,P7O6(),DtH_HANDLE,0,t1)
    
end

fig = figure('Name','Error material derivative');
plot(mw,Err(:,1),'r.-',...
    mw,Err(:,2),'go-',...
    mw,Err(:,3),'bx-',...
    mw,Err(:,4),'k+-',...
    mw,Err(:,5),'y*-',...
    mw,Err(:,6),'cs-',...
    mw,Err(:,7),'ms-'); grid('on');
set(gca,'YScale','log','XScale','log');
xlabel(['{\bf \Delta t}']);
ylabel('{\bf Error}');
legend('Location','Northwest','Lagrange','LagrangeQ','Euler+UpFEM','Interpolation','ID','DI','Euler');

[diff(log(mw)).\diff(log(Err(:,1)))...
    diff(log(mw)).\diff(log(Err(:,2)))...
    diff(log(mw)).\diff(log(Err(:,3)))...
    diff(log(mw)).\diff(log(Err(:,4)))...
    diff(log(mw)).\diff(log(Err(:,5)))...
    diff(log(mw)).\diff(log(Err(:,6)))...
    diff(log(mw)).\diff(log(Err(:,7)))]