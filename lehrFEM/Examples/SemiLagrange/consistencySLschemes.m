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

CFL=0.1;

H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)-sin(pi.*x(:,1));
D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
D12H1=@(x)-pi.*cos(pi.*x(:,1));
D22H1=@(x)0*x(:,1);

H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2).^2);
D1H1=@(x)pi*(1-x(:,2).^2).*cos(pi.*x(:,1));
D2H1=@(x)-2*x(:,2).*sin(pi.*x(:,1));
D11H1=@(x)-pi^2*(1-x(:,2).^2).*sin(pi.*x(:,1));
D12H1=@(x)-2*x(:,2).*pi.*cos(pi.*x(:,1));
D22H1=@(x)-2*sin(pi.*x(:,1));

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
    
    M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3()); 
    
    M2nd = assemMat_W1F2nd(Mesh,@MASS_W1F2nd,@(x,varargin)1, P7O6());
    W2nd = assemLoad_W1F2nd(Mesh,P7O6(),F2HANDLE);
    U2nd = assemLoad_W1F2nd(Mesh,P7O6(),H_HANDLE);
    %W = assemLoad_W1F(Mesh,P7O6(),F2HANDLE);
    %U = assemLoad_W1F(Mesh,P7O6(),F1HANDLE);
   
    % Interpolants of test and trial function
    u2nd = M2nd\U2nd;
    w2nd = M2nd\W2nd;
    w = assemCochain_1f(NewMesh,F2HANDLE,gauleg(0,1,10));
    u= assemCochain_1f(NewMesh,H_HANDLE,gauleg(0,1,10));
    
    % Mesh Width
    h(i) = get_MeshWidth(NewMesh);
    
    % SemiLagrange Quad-based  
    tau = CFL*h(i)/norm([v1,v2]);
    pbB = trace_bcenters(NewMesh, Dir_HANDLE,JvDir_Handle,-tau);
    %P_q1 = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB); 
    %P_q12nd = assemMat_SemiLagQuad_W1F_bary2nd(NewMesh, pbB);  
    P_q1 = assemMat_SemiLagQuad_W1F_strang1(NewMesh, @(x)tau*Dir_HANDLE(x),JvDir_Handle,-tau);
    P_q12nd = assemMat_SemiLagQuad_W1F_strang12nd(NewMesh, Dir_Handle,JvDir_Handle,-h);
    
    %SemiLagrange Interpol-based
    directions = Dir_HANDLE(NewMesh.Coordinates,0);
    pbV = trace_vertices(NewMesh,-tau*directions);
    P_i1 = assemMat_SemiLag_W1F(NewMesh, pbV);
    
    % SemiLagrange Patches 1st
    pbV = trace_vertices(NewMesh,tau*directions);
    NewMesh = init_LEB(NewMesh);
    NewMesh = add_Patches(NewMesh);
    NewMesh = add_Edge2Elem(NewMesh);
    defMesh = NewMesh;
    defMesh.Coordinates = pbV(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbV);
    
    P_pat = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    M_pat = assemMat_MassQuad_W1F_patches(NewMesh, defMesh,intersec);
   
    % SemiLagrange Patches 2nd
%     tau2 = CFL*h(i)^2/norm([v1,v2]); 
%     pbV = trace_vertices(NewMesh,tau2*directions);
%     NewMesh = init_LEB(NewMesh);
%     NewMesh = add_Patches(NewMesh);
%     NewMesh = add_Edge2Elem(NewMesh);
%     defMesh = NewMesh;
%     defMesh.Coordinates = pbV(:,[1 2]);
%     intersec = aff_elems2(NewMesh, defMesh,pbV);
    P_pat2nd = assemMat_SemiLagQuad_W1F_patches2nd(NewMesh, defMesh,intersec);
    M_pat2nd = assemMat_MassQuad_W1F_patches2nd(NewMesh, defMesh,intersec);
    nEdges=size(Mesh.Edges(),1);

    L1c = cDhv+cDvh+chDv;
    L1m = 0.5*cID1+cDvh+chDv;
    
    Err(i,1) = abs(w'*(M-M*P_i1)/tau*u-L1c);
    Err(i,2) = abs(w'*(M-M*P_i1)/tau*u-L1m);
    Err(i,3) = abs(w'*(M-P_q1)/tau*u-L1c);
    Err(i,4) = abs(w'*(M-P_q1)/tau*u-L1m);
    Err(i,5) = abs(w'*(M_pat-P_pat)/tau*u-L1c);
    Err(i,6) = abs(w'*(M_pat-P_pat)/tau*u-L1m);
    Err(i,7) = abs(w2nd'*(M2nd-P_q12nd)/tau*u2nd-L1c);
    Err(i,8) = abs(w2nd'*(M2nd-P_q12nd)/tau*u2nd-L1m);
    Err(i,9) = abs(w2nd'*(M_pat2nd-P_pat2nd)/tau*u2nd-L1c);
    Err(i,10) = abs(w2nd'*(M_pat2nd-P_pat2nd)/tau*u2nd-L1m)
end

fig = figure('Name','Consistency error Lie-derivatives');
plot(h,Err(:,1),'r-',...
    h,Err(:,2),'ro-',...
    h,Err(:,3),'m-',...
    h,Err(:,4),'mo-',...
    h,Err(:,5),'g-',...
    h,Err(:,6),'go-',...
    h,Err(:,7),'m--',...
    h,Err(:,8),'mo--',...
    h,Err(:,9),'g--',...
    h,Err(:,10),'go--'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf L^2-Error}');
add_Slope(gca,'NorthEast',1);
add_Slope(gca,'SouthEast',2);
legend('SL-Int-c','SL-Int-m','SL-Quad-c','SL-Quad-m','SL-Pat-c','SL-Pat-m','SL-Quad2-c','SL-Quad2-m','SL-Pat2-c','SL-Pat2-m');