% consistency of bilinarforms approximating the LieDerivative
% parts i_v d und d i_v for 1 and 2 forms, using extrusion contraction

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all
NREFS =4;  % Number of unifrom red refinements
JIG =2;

CFL=0.3;

H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)-sin(pi.*x(:,1));
D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)-pi.*cos(pi.*x(:,1));
D22H1=@(x)0*x(:,1);

% H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2).^2);
% D1H1=@(x)pi*(1-x(:,2).^2).*cos(pi.*x(:,1));
% D2H1=@(x)-2*x(:,2).*sin(pi.*x(:,1));
% D11H1=@(x)-pi^2*(1-x(:,2).^2).*sin(pi.*x(:,1));
% D12H1=@(x)-2*x(:,2).*pi.*cos(pi.*x(:,1));
% D22H1=@(x)-2*sin(pi.*x(:,1));

H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
D11H2=@(x)(1-x(:,2).^2).*(-2);
D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
D22H2=@(x)(-2).*(1-x(:,1).^2);

v1=0.66;
v2=1;
% DIR1=@(x)v1*x(2).*(1-x(:,2)).*x(1).*(1-x(:,1));
% D1DIR1=@(x)v1*x(2).*(1-x(:,2)).*(1-2*x(:,1));
% D2DIR1=@(x)v1*(1-2*x(:,2)).*x(1).*(1-x(:,1));
DIR1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

% DIR1=@(x)v1*(1-x(:,2).^2);
% D1DIR1=@(x)0*x(:,1);
% D2DIR1=@(x)v1*(-2*x(:,2));
% DIR2=@(x)0.2+sin(pi*x(:,1));
% D1DIR2=@(x)pi*cos(pi*x(:,1));
% D2DIR2=@(x)0*(x(:,2));

Dir_HANDLE=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

H_HANDLE=@(x,flag,varargin)...
    [H1(x) H2(x)];

%  -v x curl h
ID1H_HANDLE=@(x,flag,varargin)...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
    DIR1(x).*(D1H2(x)-D2H1(x))];

% v h
I1H_HANDLE=@(x,flag,varargin)...
    [DIR1(x).*H1(x)+DIR2(x).*H2(x)];

% grad v h
DI1H_HANDLE=@(x,flag,varargin)...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
    D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

% D h v
DHDIR_HANDLE=@(x,flag,varargin)...
    [D1H1(x).*DIR1(x)+D2H1(x).*DIR2(x) ...
    D1H2(x).*DIR1(x)+D2H2(x).*DIR2(x)];

% D v h
DDIRH_HANDLE=@(x,flag,varargin)...
    [D1DIR1(x).*H1(x)+D2DIR1(x).*H2(x) ...
    D1DIR2(x).*H1(x)+D2DIR2(x).*H2(x)];

%  h x curl v
HDDIR_HANDLE=@(x,flag,varargin)...
    [H2(x).*(D1DIR2(x)-D2DIR1(x)) ...
    -H1(x).*(D1DIR2(x)-D2DIR1(x))];

% v * div h
ID2H_HANDLE=@(x,flag,varargin)...
    [DIR1(x).*(D1H1(x)+D2H2(x)) ...
    DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h)
DI2H_HANDLE=@(x,flag,varargin)...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
    +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];

% exact values of bi-forms
F2HANDLE=@(x,varargin) [x(:,1).^2+x(:,2) x(:,2).*x(:,1)];
%F2HANDLE=@(x,varargin) [x(:,1).^2+x(:,2) x(:,1)];
%F2HANDLE=@(x,varargin) [ones(size(x,1),1) -x(:,1).^2+3*x(:,2) ];

F2sHANDLE=@(x,varargin) [ -x(:,1)+x(:,2).^3 ];

F1HANDLE=@(x,varargin) I1H_HANDLE(x,0);
cI1=innerproduct(F1HANDLE,F2sHANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) ID1H_HANDLE(x,0);
cID1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) DI1H_HANDLE(x,0);
cDI1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) DHDIR_HANDLE(x,0);
cDhv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) DDIRH_HANDLE(x,0);
cDvh=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) HDDIR_HANDLE(x,0);
chDv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) ID2H_HANDLE(x,0);
cID2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) DI2H_HANDLE(x,0);
cDI2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

[cID1 cDI1 cID2 cDI2 cDvh]

% Load mesh from file

% initial mesh
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');

% Add edge data structure

Mesh = add_Edges(Mesh);
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

Err=zeros(NREFS,10);
h=zeros(NREFS,1);

% refine
for i=1:NREFS
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
 
%     SP = SparsePattern(NewMesh);
%     SP=1-(SP>0);
    
    M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3());
    M_LFE = assemMat_LFE(NewMesh,@MASS_Lump_LFE);

    % Interpolants of test and trial function
    w = assemCochain_1f(NewMesh,F2HANDLE,gauleg(0,1,10));
    ws = assemCochain_0f(NewMesh,F2sHANDLE);
    u= assemCochain_1f(NewMesh,H_HANDLE,gauleg(0,1,10));
    Duv= assemCochain_1f(NewMesh,DDIRH_HANDLE,gauleg(0,1,10));

    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

    % contraction operators
    ContrOne=assemMat_Contr1f(NewMesh,Dir_HANDLE);  % contraction of one forms
    V=Dir_HANDLE(NewMesh.Coordinates);
    %ContrTwo=assemMat_Contr2f_fast(NewMesh,Dir_HANDLE);                 % contraction of two forms
    %ContrTwo=assemMat_Contr2f(NewMesh,V);                 % contraction of two forms edgewise upwind 
    ContrTwo=assemMat_Contr_1FVertUP(NewMesh,Dir_HANDLE);   % contraction of two forms( vertexwise upwind)
    % Contraction Extrusion
    
    % exterior calculus
    ID1 = M*ContrTwo*TopRot;          % (-v x curl u,w)
    %ID1=assemMat_W1F(NewMesh,@STIMA_ContrRot_UP,Dir_HANDLE);
    DI1 = M*TopGrad*ContrOne;         % (grad(v.u),w)
    DI2 = -ID1';                                 % (-curl(v x u),w)=(u,v x curl w)
    ID2 = -DI1';                                 % (v div u,w)=-(u,grad(v.w))

    % standard FEM
    %ID1_FEM = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_HANDLE, P7O6());
    % upwind Quadrature consistent
    ID1_FEM = assemMat_ContrRot_UPQuad(NewMesh,Dir_HANDLE);
    
    GC_FEM = assemMat_W1F(NewMesh,@STIMA_GradContr,Dir_HANDLE, gauleg(0,1,3));
   % LieInn_FEM = assemMat_Lie_Inn_W1F(NewMesh,Dir_HANDLE);
   % LieVol_FEM = assemMat_W1F(NewMesh,@STIMA_Lie,Dir_HANDLE,JvDir_Handle);
   % LieCompact = assemMat_W1F(NewMesh,@STIMA_LieCompact,Dir_HANDLE);

    %ContrOne_fem = assemMat_W1F(NewMesh,@STIMA_GradContrQ1, Dir_HANDLE);  % contraction of one forms
    %DI1_FEM = ContrOne_fem;
    ContrOne_fem=assemMat_Contr1f_FEM(NewMesh,Dir_HANDLE,P7O6());  % contraction of one forms
    DI1_FEM=M*TopGrad*(M_LFE\ContrOne_fem);
    DI2_FEM = -ID1_FEM';
    ID2_FEM = -DI1_FEM';
    
    % Mesh Width
    h(i) = get_MeshWidth(NewMesh);
    
    % SemiLagrange Quad-based  
    tau = CFL*h(i)/norm([v1,v2]);
    %tau=0.5*h(i);
    pbB = trace_bcenters(NewMesh, Dir_HANDLE,JvDir_Handle,-tau);
    P_q1 = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);  
    pbB = trace_bcenters(NewMesh, Dir_HANDLE,JvDir_Handle,tau);
    P_q2 = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);  
    % P = assemMat_SemiLagQuad_W1F_strang1(NewMesh, @(x)Dir_Handle(x,0),JvDir_Handle,-tau);      % strang1  

    %SemiLagrange Interpol-based
    directions = Dir_HANDLE(NewMesh.Coordinates,0);
    %pbV = trace_vertices(NewMesh,-tau*directions);
    pbV = trace_vertices_W(NewMesh,Dir_HANDLE,JvDir_Handle,-tau);
    P_i1 = assemMat_SemiLag_W1F(NewMesh, pbV);
    %pbV = trace_vertices(NewMesh,tau*directions); 
    pbV = trace_vertices_W(NewMesh,Dir_HANDLE,JvDir_Handle,tau);
    P_i2 = assemMat_SemiLag_W1F(NewMesh, pbV);

    pbV = trace_vertices(NewMesh,tau*directions);
    NewMesh = init_LEB(NewMesh);
    NewMesh = add_Patches(NewMesh);
    NewMesh = add_Edge2Elem(NewMesh);
    defMesh = NewMesh;
    defMesh.Coordinates = pbV(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbV);
    
    P_pat = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    %     P_patVol = assemMat_SemiLagQuad_W1F_patchesVol(NewMesh, defMesh,intersec);
    %     P_patInn = assemMat_SemiLagQuad_W1F_patchesInn(NewMesh, defMesh,intersec);
    M_pat = assemMat_MassQuad_W1F_patches(NewMesh, defMesh,intersec);
     
    %pbV = trace_vertices(NewMesh,-tau*directions);
    %defMesh.Coordinates = pbV(:,[1 2]);
    %ntersec = aff_elems2(NewMesh, defMesh,pbV);
    %P_pat_ad = assemMat_SemiLagQuad8_W1F_patches(NewMesh, defMesh,intersec);
    
    L1c = cDhv+cDvh+chDv;
    L1m = 0.5*cID1+cDvh+chDv;
    L2c = cDI2+cID2;
    L2m = 0.5*cDI2-cDvh;
    
    % convection discretizations
    [w'*ID2_FEM*u cID2];
    ErrC(i,1) = abs(w'*DI1*u-cDI1)/sqrt(w'*M*w);
    ErrC(i,2) = abs(w'*ID1*u-cID1)/sqrt(w'*M*w);
    ErrC(i,3) = abs(w'*DI2*u-cDI2)/sqrt(w'*M*w);
    ErrC(i,4) = abs(w'*ID2*u-cID2)/sqrt(w'*M*w);
    ErrC(i,5) = abs(w'*DI1_FEM*u-cDI1)/sqrt(w'*M*w);
    ErrC(i,6) = abs(w'*ID1_FEM*u-cID1)/sqrt(w'*M*w);
    ErrC(i,7) = abs(w'*DI2_FEM*u-cDI2)/sqrt(w'*M*w);
    ErrC(i,8) = abs(w'*ID2_FEM*u-cID2)/sqrt(w'*M*w);
%    ErrC(i,9) = abs(w'*(LieInn_FEM+LieVol_FEM)*u-L1c);
%     ErrC(i,10) =abs(w'*(LieInn_FEM+GC_FEM+ID1_FEM)*u-L1c);
%    ErrC(i,10) =abs(w'*(LieCompact)*u-L1c)
    
    Err(i,1) = abs(w'*(M-M*P_i1)/tau*u-L1c);
    Err(i,2) = abs(w'*(M-M*P_i1)/tau*u-L1m);
    Err(i,3) = abs(w'*(M-P_i2'*M)/tau*u-L2c);
    Err(i,4) = abs(w'*(M-P_i2'*M)/tau*u-L2m);
    Err(i,5) = abs(w'*(M-P_q1)/tau*u-L1c);
    Err(i,6) = abs(w'*(M-P_q1)/tau*u-L1m);
    Err(i,7) = abs(w'*(M-P_q2')/tau*u-L2c);
    Err(i,8) = abs(w'*(M-P_q2')/tau*u-L2m);
    Err(i,9) = abs(w'*(M_pat-P_pat)/tau*u-L1c);
    Err(i,10) = abs(w'*(M_pat-P_pat)/tau*u-L1m)
%  Err(i,2) = abs(w'*(M-P)/tau*u-) % is CFL-dependent
%  Err(i,3) = abs(w'*(M-P)/tau*u-cDhv-cDvh-chDv) % is CFL-dependent
end
fig = figure('Name','Consistency error convective bi-forms');
plot(h,ErrC(:,1),'r-',...
    h,ErrC(:,2),'ro-',...
    h,ErrC(:,3),'b-',...
    h,ErrC(:,4),'bo-',...
    h,ErrC(:,5),'m-',...
    h,ErrC(:,6),'mo-',...
    h,ErrC(:,7),'k-',...
    h,ErrC(:,8),'ko-',...
    h,ErrC(:,9),'g-',...
    h,ErrC(:,10),'go-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf L^2-Error}');
add_Slope(gca,'NorthEast',1);
add_Slope(gca,'SouthEast',2);
legend('DI_1','ID_1','DI_2','ID_2','DI_1-FEM','ID_1-FEM','DI_2-FEM','ID_2-FEM','L1c-FEM_p','L1m-FEM_p');

fig = figure('Name','Consistency error Lie-derivatives');
plot(h,Err(:,1),'r-',...
    h,Err(:,2),'ro-',...
    h,Err(:,3),'b-',...
    h,Err(:,4),'bo-',...
    h,Err(:,5),'m-',...
    h,Err(:,6),'mo-',...
    h,Err(:,7),'k-',...
    h,Err(:,8),'ko-',...
    h,Err(:,9),'g-',...
    h,Err(:,10),'go-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf L^2-Error}');
add_Slope(gca,'NorthEast',1);
add_Slope(gca,'SouthEast',2);
legend('SL_1-c','SL_1-m','SL_2-c','SL_2-m','SL_1FEM-c','SL_1FEM-m','SL_2FEM-c','SL_2FEM-m','SL_1FEM_p-c','SL_1FEM_p-m');

fig = figure('Name','Consistency error Lie-derivatives');
plot(h,Err(:,1),'r-',...
    h,Err(:,2),'ro-',...
    h,Err(:,5),'m-',...
    h,Err(:,6),'mo-',...
    h,Err(:,9),'g-',...
    h,Err(:,10),'go-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf L^2-Error}');
add_Slope(gca,'NorthEast',1);
add_Slope(gca,'SouthEast',2);
legend('SL_1-c','SL_1-m','SL_1FEM-c','SL_1FEM-m','SL_1FEM_p-c','SL_1FEM_p-m');


[diff(log(h)).\diff(log(ErrC(:,1)))...
    diff(log(h)).\diff(log(ErrC(:,2)))...
    diff(log(h)).\diff(log(ErrC(:,3)))...
    diff(log(h)).\diff(log(ErrC(:,4)))...
    diff(log(h)).\diff(log(ErrC(:,5)))...
    diff(log(h)).\diff(log(ErrC(:,6)))...
    diff(log(h)).\diff(log(ErrC(:,7)))...
    diff(log(h)).\diff(log(ErrC(:,8)))]

[diff(log(h)).\diff(log(Err(:,1)))...
    diff(log(h)).\diff(log(Err(:,2)))...
    diff(log(h)).\diff(log(Err(:,3)))...
    diff(log(h)).\diff(log(Err(:,4)))...
    diff(log(h)).\diff(log(Err(:,5)))...
    diff(log(h)).\diff(log(Err(:,6)))...
    diff(log(h)).\diff(log(Err(:,7)))...
    diff(log(h)).\diff(log(Err(:,8)))]