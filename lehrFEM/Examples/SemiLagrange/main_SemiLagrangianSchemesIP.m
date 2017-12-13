function main_SemiLagrangianSchemesIP(CFL,NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
filename= ['./results/example1IP_',num2str(CFL)];
%clear all;b
NREFS_init =0;     % Number of uniform red refinements
% NREFS = 3;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

%CFL =0.8; %
T0 = 0;
T1= 0.7;

mw = zeros(NREFS,1);
steps=zeros(NREFS,1);

c_CC =0;
c_Dt = 1;
c_ID = 1;
c_DI = c_ID;
    
T_Handle=@(t,varargin) 1;%cos(2*pi*t);              %1/(1+t);%1;%
DT_Handle=@(t,varargin) 0;%-2*pi*sin(2*pi*t);  %-1/(1+t)^2;%

% H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
% D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)-sin(pi.*x(:,1));
% D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)-pi.*cos(pi.*x(:,1));
% D22H1=@(x)0*x(:,1);
H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
D11H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));

H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
D11H2=@(x)(1-x(:,2).^2).*(-2);
D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
D22H2=@(x)(-2).*(1-x(:,1).^2);

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

% 
% % 2-forms: stationary solution with zero source term, div-free veolicity
B1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
D1B1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
D2B1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
D11B1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
D12B1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
D22B1=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
B2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
D1B2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
D2B2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
D11B2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
D12B2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
D22B2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);

v1=0.5;
v2=0.5;
DIR1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
D1DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
D2DIR1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
D11DIR1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
D12DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
D22DIR1=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
DIR2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
D1DIR2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
D2DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
D11DIR2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
D12DIR2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
D22DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
% 
% 1-forms: stationary solution with zero source term, div-free velocity
H1=@(x) 1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
D1H1=@(x) 1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
D2H1=@(x) -4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
D11H1=@(x) -6 * x(:,1).*(1-x(:,2).^2).^2;
D12H1=@(x) -4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
D22H1=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
H2 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
D1H2=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
D2H2=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
D11H2=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
D12H2=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
D22H2=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);

% v1=0.5;
% v2=0.5;
% DIR1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2DIR1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% D11DIR1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
% D12DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
% D22DIR1=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
% DIR2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1DIR2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
% D11DIR2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
% D12DIR2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
% D22DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
% % 
% % 1-forms: stationary solution curl-free velocity
% DIR1=@(x) 1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1DIR1=@(x) 1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2DIR1=@(x) -4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
% D11DIR1=@(x) -6 * x(:,1).*(1-x(:,2).^2).^2;
% D12DIR1=@(x) -4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
% D22DIR1=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
% DIR2 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1DIR2=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2DIR2=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% D11DIR2=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
% D12DIR2=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
% D22DIR2=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
% 
% v1=0.5;
% v2=0.5;
% H1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1H1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2H1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% D11H1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
% D12H1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
% D22H1=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
% H2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1H2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2H2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
% D11H2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
% D12H2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
% D22H2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[H1(x) H2(x)];

% - curl h
CH_Handle=@(x,flag,t,varargin)T_Handle(t).*...
    D2H1(x) +D1H2(x);

% curl curl h
CCH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22H1(x) + D12H2(x) ...
    D12H1(x)-D11H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];

%  -v x curl h
ID1H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
    DIR1(x).*(D1H2(x)-D2H1(x))];

% grad v h
DI1H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
    D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

% D h v
DHDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1H1(x).*DIR1(x)+D2H1(x).*DIR2(x) ...
    D1H2(x).*DIR1(x)+D2H2(x).*DIR2(x)];

% D v h
DDIRH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+D2DIR1(x).*H2(x) ...
    D1DIR2(x).*H1(x)+D2DIR2(x).*H2(x)];

% div v h
DivDIRH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [(D1DIR1(x)+D2DIR2(x)).*H1(x) ...
    (D1DIR1(x)+D2DIR2(x)).*H2(x)];

%  h x curl v
HDDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [H2(x).*(D1DIR2(x)-D2DIR1(x)) ...
    -H1(x).*(D1DIR2(x)-D2DIR1(x))];

% v * div h
ID2H_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [DIR1(x).*(D1H1(x)+D2H2(x)) ...
    DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h)
DI2H_HANDLE=@(x,flag,t,vvarargin)T_Handle(t).*...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
     +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];

F_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    c_ID * ID1H_HANDLE(x,flag,t) + ...
    c_DI * DI1H_HANDLE(x,flag,t);

% F_HANDLE=@(x,flag,t,varargin)...
%     c_CC * CCH_HANDLE(x,flag,t) + ...
%     c_Dt * DtH_HANDLE(x,flag,t) + ...
%     HDDIR_HANDLE(x,flag,t) + ...
%     DDIRH_HANDLE(x,flag,t) + ...
%     DHDIR_HANDLE(x,flag,t);

B_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[B1(x) B2(x)];

% - curl h
CB_Handle=@(x,flag,t,varargin)T_Handle(t).*...
    D2B1(x) +D1B2(x);

% curl curl h
CCB_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22B1(x) + D12B2(x) ...
    D12B1(x)-D11B2(x)];

% d_t h
DtB_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ B1(x) B2(x)];

%  -v x curl h
ID1B_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-DIR2(x).*(D1B2(x)-D2B1(x)) ...
    DIR1(x).*(D1B2(x)-D2B1(x))];

% grad v h
DI1B_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*B1(x)+DIR1(x).*D1B1(x)+D1DIR2(x).*B2(x)+DIR2(x).*D1B2(x) ...
    D2DIR1(x).*B1(x)+DIR1(x).*D2B1(x)+D2DIR2(x).*B2(x)+DIR2(x).*D2B2(x)];

% D h v
DBDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1B1(x).*DIR1(x)+D2B1(x).*DIR2(x) ...
    D1B2(x).*DIR1(x)+D2B2(x).*DIR2(x)];

% D v h
DDIRB_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*B1(x)+D2DIR1(x).*B2(x) ...
    D1DIR2(x).*B1(x)+D2DIR2(x).*B2(x)];

% div v h
DivDIRH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [(D1DIR1(x)+D2DIR2(x)).*B1(x) ...
    (D1DIR1(x)+D2DIR2(x)).*B2(x)];

%  h x curl v
HDDIR_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [B2(x).*(D1DIR2(x)-D2DIR1(x)) ...
    -B1(x).*(D1DIR2(x)-D2DIR1(x))];

% v * div h
ID2B_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [DIR1(x).*(D1B1(x)+D2B2(x)) ...
    DIR2(x).*(D1B1(x)+D2B2(x))];

% -curl(v x h)
DI2B_HANDLE=@(x,flag,t,vvarargin)T_Handle(t).*...
    [-D2DIR1(x).*B2(x)-DIR1(x).*D2B2(x)+D2DIR2(x).*B1(x)+DIR2(x).*D2B1(x) ...
     +D1DIR1(x).*B2(x)+DIR1(x).*D1B2(x)-D1DIR2(x).*B1(x)-DIR2(x).*D1B1(x)];

G_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCB_HANDLE(x,flag,t) + ...
    c_Dt * DtB_HANDLE(x,flag,t) + ...
    ID2B_HANDLE(x,flag,t) + ...
    DI2B_HANDLE(x,flag,t);

% G_HANDLE=@(x,flag,t,varargin)...
%     c_CC * CCH_HANDLE(x,flag,t) + ...
%     c_Dt * DtH_HANDLE(x,flag,t) - ...
%     DDIRH_HANDLE(x,flag,t) + ...
%     0.5*DI2H_HANDLE(x,flag,t);

% Load mesh from file
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

% refine
for i=1:NREFS_init
    Mesh=refine_REG(Mesh);
end

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

    mw(j) = get_MeshWidth(NewMesh);

    h = CFL*mw(j)/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps;
    steps(j) = nsteps;

    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
    C = assemMat_W1F(NewMesh,@STIMA_Curl_W1F,SIGMA_HANDLE,P7O6());
    %M = assemMat_Mass1fD(NewMesh);
    
    % discrete pullback
    % Interpolation
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
    %pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h);
    %          direction = Dir_Handle(NewMesh.Coordinates,0);
    %          pbVp = trace_vertices(NewMesh,h*direction);
    P_i_dir = assemMat_SemiLag_W1F(NewMesh, pbVm);
    %P_i_ad = assemMat_SemiLag_W1F(NewMesh, pbVp);
    
    % transported Mesh
    pMesh.Coordinates=pbVm(:,[1,2]);
    pMesh.Elements=Mesh.Elements;
    pMesh = add_Edges(pMesh);
    Loc = get_BdEdges(Mesh);
    pMesh.BdFlags = zeros(size(pMesh.Edges,1),1);
    pMesh.BdFlags(Loc) = -1;
    Mm = assemMat_W1F(pMesh,@MASS_W1F,MU_HANDLE, P3O3());
    
%     AvJump=assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_AvgJumpW1F,Dir_Handle,gauleg(0,1,5));
%     
%     Lie1 = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle,P7O6());
%     Lie2 =  assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,Dir_Handle,gauleg(0,1,10));
%     Lie2bnd = assemMat_Bnd_LieW1F(NewMesh,[-2],@STIMA_Bnd_LieW1F,Dir_Handle,gauleg(0,1,5));
%     Lie2pen = assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,Dir_Handle,gauleg(0,1,10),0.5);
%     Lie = Lie1 + Lie2+Lie2bnd-Lie2pen;
  
    % timestepping

    %[interpolation quadrature patches]
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old = assemCochain_1f(NewMesh,B_Handle,gauleg(0,1,10),0);
    plot_Norm_W1F(B_old,NewMesh);colorbar;
    
    time(1,j) = 0;
    L2ErrH(1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,2) = L2Err_W1F(NewMesh,B_old(:,1),P7O6(),B_Handle,0,0);
    
%     defMesh = NewMesh;
%     defMesh.Coordinates = pbVm(:,[1 2]);
%     intersec = aff_elems2(NewMesh, defMesh,pbVm);
%     P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);

    for i = 1:nsteps
        [i nsteps]
        time(i+1,j) = h+time(i,j);

        %Dirichlet data
        [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),(i)*h);
        H_new=[H_Dir];
        %system matrix
        A = M; 
        %righthand side
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),i*h);
        L=[Lf];
        L = h*M*L+M*P_i_dir* H_old-A*H_new;
        %solve system
        H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
        %Update vectors
        H_old=H_new;
        L2ErrH(i+1,j,1) = L2Err_W1F(NewMesh,H_old,P7O6(),H_Handle,0,i*h);
        
        %Dirichlet data
        [B_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,B_Handle,gauleg(0,1,5),(i)*h);
        B_new=[B_Dir];
        %system matrix
        A = M;
        %righthand side
        Lg = assemCochain_1f(NewMesh, G_HANDLE,gauleg(0,1,10),i*h);
        L=Lg;
        %plot_Norm_W1F(Lg,NewMesh); colorbar;
        L = h*M*L+Mm*P_i_dir*B_old-A*B_new;
        %L = P_i_ad'*M* B_old-A*B_new;
        %solve system
        %B_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
        B_new = A\L;
        
        %Update vectors
        B_old=B_new;
        L2ErrH(i+1,j,2) = L2Err_W1F(NewMesh,B_old,P7O6(),B_Handle,0,i*h);
        %plot_Norm_W1F(B_old,NewMesh);colorbar;
        H_ref = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),i*h);
        B_ref = assemCochain_1f(NewMesh,B_Handle,gauleg(0,1,10),i*h);
   
    
    end % timestep
    plot_Norm_W1F(H_old,NewMesh);colorbar;
    plot_Norm_W1F(H_ref,NewMesh);colorbar;
    
    plot_Norm_W1F(B_old,NewMesh);colorbar;
    plot_Norm_W1F(B_ref,NewMesh);colorbar;
    
    
end%  Mesh refinement

for j=1:NREFS
   j
   stopTimeL2ErrH(j,:)=L2ErrH(steps(j)+1,j,:)
end
L2ErrH(:,:,2)

figure;
hold on;
plot(mw,stopTimeL2ErrH(:,1),'r-',mw,stopTimeL2ErrH(:,2),'g-','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('direct','adjoint','Location','Northwest');
saveas(gcf,[filename,'error.fig'])
print('-depsc' ,[filename,'error.eps'])

slope=[diff(log(stopTimeL2ErrH(:,1)))./diff(log(mw)),diff(log(stopTimeL2ErrH(:,2)))./diff(log(mw))]

save([filename,'.mat'],'slope','mw','stopTimeL2ErrH','CFL','H1','H2','DIR1','DIR2')
save([filename,'.txt']','slope','mw','stopTimeL2ErrH','CFL','-ASCII')
clear all
