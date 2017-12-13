%function [time,mass,lm]=SemiLagrangeDirect(CFL)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;
clear Mesh;
%clear all;
NREFS_init = 0;     % Number of uniform red refinements
NREFS = 5;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

CFL =0.8; %
T0 = 0;
T1= 0.4;

mw = zeros(NREFS,1);
Err = [];
time = [];
div=[];
mass=[];
steps=zeros(NREFS,1);

c_CC =0;
c_Dt = 1;
c_ID = 1;
c_DI = c_ID;

T_Handle=@(t,varargin) cos(2*pi*t);             %1/(1+t);%1;%
DT_Handle=@(t,varargin) -2*pi*sin(2*pi*t);  %-1/(1+t)^2;%

H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)-sin(pi.*x(:,1));
D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)-pi.*cos(pi.*x(:,1));
D22H1=@(x)0*x(:,1);

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
DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)v2*(1-x(:,1));
% D1DIR2=@(x)-v2*ones(size(x,1),1);
% D2DIR2=@(x)0*(x(:,2));

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


% F_HANDLE=@(x,flag,t,varargin)...
%     c_CC * CCH_HANDLE(x,flag,t) + ...
%     c_Dt * DtH_HANDLE(x,flag,t) + ...
%     c_ID * ID1H_HANDLE(x,flag,t) + ...
%     c_DI * DI1H_HANDLE(x,flag,t);

F_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    HDDIR_HANDLE(x,flag,t) + ...
    DDIRH_HANDLE(x,flag,t) + ...
    DHDIR_HANDLE(x,flag,t);

F_mod_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    HDDIR_HANDLE(x,flag,t) + ...
    DDIRH_HANDLE(x,flag,t) + ...
    0.5*ID1H_HANDLE(x,flag,t);

% F_modL_HANDLE=@(x,flag,t,varargin)...
%     c_CC * CCH_HANDLE(x,flag,t) + ...
%     c_Dt * DtH_HANDLE(x,flag,t) + ...
%     ID1H_HANDLE(x,flag,t);

G_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    ID2H_HANDLE(x,flag,t) + ...
    DI2H_HANDLE(x,flag,t);

G_mod_HANDLE=@(x,flag,t,varargin)...
    c_CC * CCH_HANDLE(x,flag,t) + ...
    c_Dt * DtH_HANDLE(x,flag,t) + ...
    0.5*DI2H_HANDLE(x,flag,t) - ...
    DDIRH_HANDLE(x,flag,t);

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

    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % Assemble Curl-curl matrix, MASS matrix and load vector
    C = assemMat_W1F(NewMesh,@STIMA_Curl_W1F,SIGMA_HANDLE,P7O6());
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());

    % topological derivatives
    %TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    %TopRot=assemMat_TopRot(NewMesh);     % topological Rotation

    % discrete pullback
    % Interpolation
    directions = Dir_Handle(NewMesh.Coordinates,i*h);
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
    %pbV = trace_vertices(NewMesh,-h*directions);
    P_i = assemMat_SemiLag_W1F(NewMesh, pbVm);
    pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h);
    %pbV = trace_vertices(NewMesh,h*directions);
    P_i_ad = assemMat_SemiLag_W1F(NewMesh, pbVp);

    % Quadrature
    pbB = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
    P_qbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);
    pbB = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,h);
    P_q_adbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);

    %plot_transportMesh(NewMesh,pbV,'tas')
    %pbV = trace_vertices(NewMesh,h*directions);
    NewMesh = init_LEB(NewMesh);
    NewMesh = add_Patches(NewMesh);
    NewMesh = add_Edge2Elem(NewMesh);
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    
    P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    %pbV = trace_vertices(NewMesh,-h*directions);
    defMesh.Coordinates = pbVm(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVm);
    P_q_ad = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    % test difference of P_i and P_q for small timesteps
 
    %Velocity tensor
    %     VT = assemMat_W1F(NewMesh,@STIMA_VTensor,JvDir_Handle);

    % timestepping

    %[interpolation quadrature interpolationMod quadratureMod]
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    H_old=[H_old H_old H_old H_old H_old H_old];

    %     HL_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    %     HL_old=[H_old H_old];

    B_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old=[B_old B_old B_old B_old B_old B_old];

    time(1,j) = 0;
    L2ErrH(1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,0);
    L2ErrH(1,j,6) = L2Err_W1F(NewMesh,H_old(:,6),P7O6(),H_Handle,0,0);
    HCurlErrH(1,j,1) = HCurlSErr_W1F(NewMesh,H_old(:,1),P7O6(),CH_Handle,0,0);
    HCurlErrH(1,j,2) = HCurlSErr_W1F(NewMesh,H_old(:,2),P7O6(),CH_Handle,0,0);
    HCurlErrH(1,j,3) = HCurlSErr_W1F(NewMesh,H_old(:,3),P7O6(),CH_Handle,0,0);
    HCurlErrH(1,j,4) = HCurlSErr_W1F(NewMesh,H_old(:,4),P7O6(),CH_Handle,0,0);
    HCurlErrH(1,j,5) = HCurlSErr_W1F(NewMesh,H_old(:,5),P7O6(),CH_Handle,0,0);
    HCurlErrH(1,j,6) = HCurlSErr_W1F(NewMesh,H_old(:,6),P7O6(),CH_Handle,0,0);

    %     L2ErrHL(1,j,1) = L2Err_W1F(NewMesh,HL_old(:,1),P7O6(),H_Handle,0,0);
    %     L2ErrHL(1,j,2) = L2Err_W1F(NewMesh,HL_old(:,2),P7O6(),H_Handle,0,0);
    %     HCurlErrHL(1,j,1) = HCurlSErr_W1F(NewMesh,HL_old(:,1),P7O6(),CH_Handle,0,0);
    %     HCurlErrHL(1,j,2) = HCurlSErr_W1F(NewMesh,HL_old(:,2),P7O6(),CH_Handle,0,0);

    L2ErrB(1,j,1) = L2Err_W1F(NewMesh,B_old(:,1),P7O6(),H_Handle,0,0);
    L2ErrB(1,j,2) = L2Err_W1F(NewMesh,B_old(:,2),P7O6(),H_Handle,0,0);
    L2ErrB(1,j,3) = L2Err_W1F(NewMesh,B_old(:,3),P7O6(),H_Handle,0,0);
    L2ErrB(1,j,4) = L2Err_W1F(NewMesh,B_old(:,4),P7O6(),H_Handle,0,0);
    L2ErrB(1,j,5) = L2Err_W1F(NewMesh,B_old(:,5),P7O6(),H_Handle,0,0);
    L2ErrB(1,j,6) = L2Err_W1F(NewMesh,B_old(:,6),P7O6(),H_Handle,0,0);
    HCurlErrB(1,j,1) = HCurlSErr_W1F(NewMesh,B_old(:,1),P7O6(),CH_Handle,0,0);
    HCurlErrB(1,j,2) = HCurlSErr_W1F(NewMesh,B_old(:,2),P7O6(),CH_Handle,0,0);
    HCurlErrB(1,j,3) = HCurlSErr_W1F(NewMesh,B_old(:,3),P7O6(),CH_Handle,0,0);
    HCurlErrB(1,j,4) = HCurlSErr_W1F(NewMesh,B_old(:,4),P7O6(),CH_Handle,0,0);
    HCurlErrB(1,j,5) = HCurlSErr_W1F(NewMesh,B_old(:,5),P7O6(),CH_Handle,0,0);
    HCurlErrB(1,j,6) = HCurlSErr_W1F(NewMesh,B_old(:,6),P7O6(),CH_Handle,0,0);
    
    %  test_Patches(NewMesh,defMesh,intersec);

    for i = 1:nsteps
        [i nsteps]
        time(i+1,j) = h+time(i,j);

        %Dirichlet data
        [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),(i)*h);
        H_new=[H_Dir H_Dir H_Dir H_Dir H_Dir H_Dir];

                %[HL_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),(i)*h);
                %HL_new=[HL_Dir HL_Dir];

        [B_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),(i)*h);
        B_new=[B_Dir B_Dir B_Dir B_Dir B_Dir B_Dir];

        %system matrix
        A = c_Dt*M + c_CC*h*C;
                %AT_exp = 2*c_Dt*M + c_CC*h*C;
                %AT_imp =2*c_Dt*M + c_CC*h*C-2*h*VT;

        %righthand side
        Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),(i)*h);
        Lm = assemCochain_1f(NewMesh, F_mod_HANDLE,gauleg(0,1,10),(i)*h);
        L=[Lf Lf Lf Lm Lm Lm];
        L = h*M*L+[M*P_i*H_old(:,1) P_qbary*H_old(:,2) P_q*H_old(:,3) M*P_i*H_old(:,4) P_qbary*H_old(:,5) P_q*H_old(:,6)]-A*H_new;
        
        %LL = assemCochain_1f(NewMesh, F_modL_HANDLE,gauleg(0,1,10),(i)*h);
        %LL_exp = h*M*LL+2*P_q*HL_old(:,1)+2*h*VT*HL_old(:,1) -AT_exp*HL_new(:,1);
        %LL_imp = h*M*LL+2*P_q*HL_old(:,2)-AT_imp*HL_new(:,2);

        Gf = assemCochain_1f(NewMesh, G_HANDLE,gauleg(0,1,10),(i)*h);
        Gm = assemCochain_1f(NewMesh, G_mod_HANDLE,gauleg(0,1,10),(i)*h);
        G=[Gf Gf Gf Gm Gm Gm];
        G = h*M*G+[P_i_ad'*M*B_old(:,1) P_q_adbary'*B_old(:,2) P_q_ad'*B_old(:,3) P_i_ad'*M*B_old(:,4) P_q_adbary'*B_old(:,5) P_q_ad'*B_old(:,6)]-A*B_new;

        %solve system
        H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
        %        HL_new(FreeDofs,1) = AT_exp(FreeDofs,FreeDofs)\LL_exp(FreeDofs);
        %        HL_new(FreeDofs,2) = AT_imp(FreeDofs,FreeDofs)\LL_imp(FreeDofs);
        B_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\G(FreeDofs,:);

        %Update vectors
        H_old=H_new;
        %        HL_old=HL_new;
        B_old=B_new;

        L2ErrH(i+1,j,1) = L2Err_W1F(NewMesh,H_old(:,1),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,2) = L2Err_W1F(NewMesh,H_old(:,2),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,3) = L2Err_W1F(NewMesh,H_old(:,3),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,4) = L2Err_W1F(NewMesh,H_old(:,4),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,5) = L2Err_W1F(NewMesh,H_old(:,5),P7O6(),H_Handle,0,i*h);
        L2ErrH(i+1,j,6) = L2Err_W1F(NewMesh,H_old(:,6),P7O6(),H_Handle,0,i*h);
        HCurlErrH(i+1,j,1) = HCurlSErr_W1F(NewMesh,H_old(:,1),P7O6(),CH_Handle,0,i*h);
        HCurlErrH(i+1,j,2) = HCurlSErr_W1F(NewMesh,H_old(:,2),P7O6(),CH_Handle,0,i*h);
        HCurlErrH(i+1,j,3) = HCurlSErr_W1F(NewMesh,H_old(:,3),P7O6(),CH_Handle,0,i*h);
        HCurlErrH(i+1,j,4) = HCurlSErr_W1F(NewMesh,H_old(:,4),P7O6(),CH_Handle,0,i*h);
        HCurlErrH(i+1,j,5) = HCurlSErr_W1F(NewMesh,H_old(:,5),P7O6(),CH_Handle,0,i*h);
        HCurlErrH(i+1,j,6) = HCurlSErr_W1F(NewMesh,H_old(:,6),P7O6(),CH_Handle,0,i*h);
        
        %L2ErrHL(i+1,j,1) = L2Err_W1F(NewMesh,HL_old(:,1),P7O6(),H_Handle,0,i*h);
        %L2ErrHL(i+1,j,2) = L2Err_W1F(NewMesh,HL_old(:,2),P7O6(),H_Handle,0,i*h);
        %HCurlErrHL(i+1,j,1) = HCurlSErr_W1F(NewMesh,HL_old(:,1),P7O6(),CH_Handle,0,i*h);
        %HCurlErrHL(i+1,j,2) = HCurlSErr_W1F(NewMesh,HL_old(:,2),P7O6(),CH_Handle,0,i*h);
        
        L2ErrB(i+1,j,1) = L2Err_W1F(NewMesh,B_old(:,1),P7O6(),H_Handle,0,i*h);
        L2ErrB(i+1,j,2) = L2Err_W1F(NewMesh,B_old(:,2),P7O6(),H_Handle,0,i*h);
        L2ErrB(i+1,j,3) = L2Err_W1F(NewMesh,B_old(:,3),P7O6(),H_Handle,0,i*h);
        L2ErrB(i+1,j,4) = L2Err_W1F(NewMesh,B_old(:,4),P7O6(),H_Handle,0,i*h);
        L2ErrB(i+1,j,5) = L2Err_W1F(NewMesh,B_old(:,5),P7O6(),H_Handle,0,i*h);
        L2ErrB(i+1,j,6) = L2Err_W1F(NewMesh,B_old(:,6),P7O6(),H_Handle,0,i*h);
        HCurlErrB(i+1,j,1) = HCurlSErr_W1F(NewMesh,B_old(:,1),P7O6(),CH_Handle,0,i*h);
        HCurlErrB(i+1,j,2) = HCurlSErr_W1F(NewMesh,B_old(:,2),P7O6(),CH_Handle,0,i*h);
        HCurlErrB(i+1,j,3) = HCurlSErr_W1F(NewMesh,B_old(:,3),P7O6(),CH_Handle,0,i*h);
        HCurlErrB(i+1,j,4) = HCurlSErr_W1F(NewMesh,B_old(:,4),P7O6(),CH_Handle,0,i*h);
        HCurlErrB(i+1,j,5) = HCurlSErr_W1F(NewMesh,B_old(:,5),P7O6(),CH_Handle,0,i*h);
        HCurlErrB(i+1,j,6) = HCurlSErr_W1F(NewMesh,B_old(:,6),P7O6(),CH_Handle,0,i*h);

   end % timestep

end%  Mesh refinement

for j=1:NREFS
   stopTimeL2ErrH(j,:)=L2ErrH(steps(j)+1,j,:);
   stopTimeHCurlErrH(j,:)=HCurlErrH(steps(j)+1,j,:);
   % stopTimeL2ErrHL(j,:)=L2ErrHL(steps(j)+1,j,:);
   % stopTimeHCurlErrHL(j,:)=HCurlErrHL(steps(j)+1,j,:);
   stopTimeL2ErrB(j,:)=L2ErrB(steps(j)+1,j,:);
   stopTimeHCurlErrB(j,:)=HCurlErrB(steps(j)+1,j,:);
end

figure;

subplot(2,2,1);
%figure;
title('Direct formulation, complete Lie-derivative');
hold on;
plot(mw,stopTimeL2ErrH(:,1),'rx-',mw,stopTimeL2ErrH(:,2),'bx-',mw,stopTimeL2ErrH(:,3),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf L^2-error}');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,1)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,2)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,3)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');
%saveas(gcf,'f.fig')
%print -depsc f.eps

subplot(2,2,2);
%figure;
title('Direct formulation, incomplete Lie-derivative')
hold on;
plot(mw,stopTimeL2ErrH(:,4),'rx-',mw,stopTimeL2ErrH(:,5),'bx-',mw,stopTimeL2ErrH(:,6),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf L^2-error}');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,4)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,5)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeL2ErrH(:,6)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');
%saveas(gcf,'g.fig')
%print -depsc g.eps

subplot(2,2,3);
%figure;
title('Adjoint formulation, complete Lie-derivative')
hold on;
plot(mw,stopTimeL2ErrB(:,1),'rx-',mw,stopTimeL2ErrB(:,2),'bx-',mw,stopTimeL2ErrB(:,3),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf L^2-error}');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,1)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,2)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,3)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');

subplot(2,2,4);
%figure;
title('Adjoint formulation, incomplete Lie-derivative')
hold on;
plot(mw,stopTimeL2ErrB(:,4),'rx-',mw,stopTimeL2ErrB(:,5),'bx-',mw,stopTimeL2ErrB(:,6),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf L^2-error}');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,4)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,5)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeL2ErrB(:,6)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');

figure;

subplot(2,2,1);
%figure;
title('Direct formulation, complete Lie-derivative')
hold on;
plot(mw,stopTimeHCurlErrH(:,1),'rx-',mw,stopTimeHCurlErrH(:,2),'bx-',mw,stopTimeHCurlErrH(:,3),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf H(curl)-semi error}');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,1)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,2)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,3)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');
%saveas(gcf,'f.fig')
%print -depsc f.eps

subplot(2,2,2);
%figure;
title('Direct formulation, incomplete Lie-derivative')
hold on;
plot(mw,stopTimeHCurlErrH(:,4),'rx-',mw,stopTimeHCurlErrH(:,5),'bx-',mw,stopTimeHCurlErrH(:,6),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf H(curl)-semi error}');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,4)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,5)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeHCurlErrH(:,6)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');
%saveas(gcf,'f.fig')
%print -depsc f.eps

subplot(2,2,3);
%figure;
title('Adjoint formulation, complete Lie-derivative')
hold on;
plot(mw,stopTimeHCurlErrB(:,1),'rx-',mw,stopTimeHCurlErrB(:,2),'bx-',mw,stopTimeHCurlErrB(:,3),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf H(curl)-semi error }');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,1)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,2)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,3)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');

subplot(2,2,4);
%figure;
title('Adjoint formulation, complete Lie-derivative')
hold on;
plot(mw,stopTimeHCurlErrB(:,4),'rx-',mw,stopTimeHCurlErrB(:,5),'bx-',mw,stopTimeHCurlErrB(:,6),'gx-');
grid('on');
set(gca,'YScale','log','XScale','log');
xlabel('{\bf mesh width}');
ylabel('{\bf H(curl)-semi error }');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,4)),1);
add_Slope(gca,'East',p(1),'r');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,5)),1);
add_Slope(gca,'SouthEast',p(1),'b');
p = polyfit(log(mw),log(stopTimeHCurlErrB(:,6)),1);
add_Slope(gca,'NorthEast',p(1),'g');
legend('Interpolation','Barycenters','Patches','Location','Northwest');

% figure;
% title('direct formulation, -v x curl u')
% hold on;
% plot(mw,stopTimeL2ErrHL(:,1),'rx-',mw,stopTimeL2ErrHL(:,2),'bx-');
% grid('on');
% set(gca,'YScale','log','XScale','log');
% xlabel('{\bf mesh width}');
% ylabel('{\bf L^2-error}');
% p = polyfit(log(mw),log(stopTimeL2ErrHL(:,1)),1);
% add_Slope(gca,'East',p(1),'r');
% p = polyfit(log(mw),log(stopTimeL2ErrHL(:,2)),1);
% add_Slope(gca,'SouthEast',p(1),'b');
% legend('explizit','implizit');
%
% figure;
% title('direct formulation, -v x curl u')
% hold on;
% plot(mw,stopTimeHCurlErrHL(:,1),'rx-',mw,stopTimeHCurlErrHL(:,2),'bx-');
% grid('on');
% set(gca,'YScale','log','XScale','log');
% xlabel('{\bf mesh width}');
% ylabel('{\bf H(curl)-semi -error}');
% p = polyfit(log(mw),log(stopTimeHCurlErrHL(:,1)),1);
% add_Slope(gca,'East',p(1),'r');
% p = polyfit(log(mw),log(stopTimeHCurlErrHL(:,2)),1);
% add_Slope(gca,'SouthEast',p(1),'b');
% legend('explizit','implizit');


%saveas(gcf,'f.fig')
%print -depsc f.eps
%saveas(gcf,'g.fig')
%print -depsc g.eps
%
% fig = figure('Name','Error material derivative, Interpol');
% leg=[];
% for j=1:NREFS
%     plot(time(1:steps(j)+1,j),L2ErrH(1:steps(j)+1,j,1));
%     hold on;
%     lm(j,:)=['h=',sprintf('%1.3f',mw(j)),...
%         ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps(j)),...
%         ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw(j)*steps(j)))];
% end
% grid('on');
% set(gca,'YScale','log');
% Markers = '.ox+*sdv^<>ph';
% noMarkers = length(Markers);
% Colors = 'gbrcmy';
% noColors = length(Colors);
% H1c = findobj(gca,'Type','line');
% linehandles = [H1c];
% for K = 1:length(linehandles)
%     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
%     set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
% end
% legend(lm);
% xlabel('{\bf time}');
% ylabel('{\bf error}');
% hold off;
%
% fig = figure('Name','Error material derivative, Quad');
% leg=[];
% for j=1:NREFS
%     plot(time(1:steps(j)+1,j),L2ErrH(1:steps(j)+1,j,2));
%     hold on;
%     lm(j,:)=['h=',sprintf('%1.3f',mw(j)),...
%         ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps(j)),...
%         ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw(j)*steps(j)))];
% end
% grid('on');
% set(gca,'YScale','log');
% Markers = '.ox+*sdv^<>ph';
% noMarkers = length(Markers);
% Colors = 'gbrcmy';
% noColors = length(Colors);
% H1c = findobj(gca,'Type','line');
% linehandles = [H1c];
% for K = 1:length(linehandles)
%     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
%     set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
% end
% legend(lm);
% xlabel('{\bf time}');
% ylabel('{\bf error}');
% hold off;
%
% fig = figure('Name','Error material derivative, Interpol');
% leg=[];
% for j=1:NREFS
%     plot(time(1:steps(j)+1,j),L2ErrH(1:steps(j)+1,j,3));
%     hold on;
%     lm(j,:)=['h=',sprintf('%1.3f',mw(j)),...
%         ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps(j)),...
%         ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw(j)*steps(j)))];
% end
% grid('on');
% set(gca,'YScale','log');
% Markers = '.ox+*sdv^<>ph';
% noMarkers = length(Markers);
% Colors = 'gbrcmy';
% noColors = length(Colors);
% H1c = findobj(gca,'Type','line');
% linehandles = [H1c];
% for K = 1:length(linehandles)
%     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
%     set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
% end
% legend(lm);
% xlabel('{\bf time}');
% ylabel('{\bf error}');
% hold off;
%
% fig = figure('Name','Error material derivative, Quad');
% leg=[];
% for j=1:NREFS
%     plot(time(1:steps(j)+1,j),L2ErrH(1:steps(j)+1,j,4));
%     hold on;
%     lm(j,:)=['h=',sprintf('%1.3f',mw(j)),...
%         ' \Delta t=',sprintf('%1.3f',(T1-T0)/steps(j)),...
%         ' CFL=',sprintf('%1.3f',norm([v1 v2])*(T1-T0)/(mw(j)*steps(j)))];
% end
% grid('on');
% set(gca,'YScale','log');
% Markers = '.ox+*sdv^<>ph';
% noMarkers = length(Markers);
% Colors = 'gbrcmy';
% noColors = length(Colors);
% H1c = findobj(gca,'Type','line');
% linehandles = [H1c];
% for K = 1:length(linehandles)
%     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
%     set(linehandles(K),'Color',Colors(1+mod(K-1,noColors)));
% end
% legend(lm);
% xlabel('{\bf time}');
% ylabel('{\bf error}');
% hold off;
