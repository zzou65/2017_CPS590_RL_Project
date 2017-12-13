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
JIG =1;

CFL=0.7;

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
DIR1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));


Dir_HANDLE=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];



H_HANDLE=@(x,flag,varargin)...
    [H1(x) H2(x)];

EH_HANDLE=@(x,flag,h,varargin)...
    [H1(x) H2(x)] - ...
    [H1(x-h*Dir_HANDLE(x,0)) H2(x-h*Dir_HANDLE(x,0))] - ...
    h*[D1DIR1(x).*H1(x-h*Dir_HANDLE(x,0))+D1DIR2(x).*H2(x-h*Dir_HANDLE(x,0))...
     D2DIR1(x).*H1(x-h*Dir_HANDLE(x,0))+D2DIR2(x).*H2(x-h*Dir_HANDLE(x,0))];

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
%F2HANDLE=@(x,varargin) [x(:,1).^2+x(:,2) x(:,2).*x(:,1)];
%F2HANDLE=@(x,varargin) [x(:,1).^2+x(:,2) x(:,1)];
%F2HANDLE=@(x,varargin) [ones(size(x,1),1) -x(:,1).^2+3*x(:,2) ];

% F2sHANDLE=@(x,varargin) [ -x(:,1)+x(:,2).^3 ];
% 
% F1HANDLE=@(x,varargin) I1H_HANDLE(x,0);
% cI1=innerproduct(F1HANDLE,F2sHANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) ID1H_HANDLE(x,0);
% cID1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) DI1H_HANDLE(x,0);
% cDI1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) DHDIR_HANDLE(x,0);
% cDhv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) DDIRH_HANDLE(x,0);
% cDvh=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) HDDIR_HANDLE(x,0);
% chDv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) ID2H_HANDLE(x,0);
% cID2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) DI2H_HANDLE(x,0);
% cDI2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% [cID1 cDI1 cID2 cDI2 cDvh]

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
  
    % Interpolants of trial function
    u= assemCochain_1f(NewMesh,H_HANDLE,gauleg(0,1,10));
    
    % Mesh Width
    h(i) = get_MeshWidth(NewMesh);
    
    % SemiLagrange Quad-based  
    tau = 10^-3 ;%CFL*h(i)/norm([v1,v2]);%
    
    %SemiLagrange Interpol-based
    directions = Dir_HANDLE(NewMesh.Coordinates,0);
    %pbV = trace_vertices(NewMesh,-tau*directions);
    pbV = trace_vertices_W(NewMesh,Dir_HANDLE,JvDir_Handle,-tau);
    P_i1 = assemMat_SemiLag_W1F(NewMesh, pbV);

    pbV = trace_vertices(NewMesh,tau*directions);
    NewMesh = init_LEB(NewMesh);
    NewMesh = add_Patches(NewMesh);
    NewMesh = add_Edge2Elem(NewMesh);
    defMesh = NewMesh;
    defMesh.Coordinates = pbV(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbV);
    
    P_pat = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    M_pat = assemMat_MassQuad_W1F_patches(NewMesh, defMesh,intersec);
    M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3());
    
    Err(i,1) = L2Err_W1F(NewMesh,u-P_i1*u,P7O6(),EH_HANDLE,tau,0);
    Err(i,2) = L2Err_W1F(NewMesh,u-M\(P_pat*u),P7O6(),EH_HANDLE,tau,0);
end

fig = figure('Name','Consistency error Lie-derivatives');
plot(h,Err(:,1),'r-',...
    h,Err(:,2),'bo-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf L^2-Error}');
add_Slope(gca,'NorthEast',1);
add_Slope(gca,'SouthEast',2);
legend('Interpol','Patches');