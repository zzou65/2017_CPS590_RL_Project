% consistency of bilinarforms approximating the LieDerivative
% parts i_v d und d i_v for 1 and 2 forms, using extrusion contraction

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all
NREFS =6;  % Number of unifrom red refinements
JIG =2;

%CFL=0.1;

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
% DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,1));
D1DIR2=@(x)-v2*ones(size(x,1),1);
D2DIR2=@(x)0*(x(:,2));

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

% F1HANDLE=@(x,varargin) I1H_HANDLE(x,0);
% cI1=innerproduct(F1HANDLE,F2sHANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) ID1H_HANDLE(x,0);
cID1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

 F1HANDLE=@(x,varargin) DI1H_HANDLE(x,0);
cDI1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
F1HANDLE=@(x,varargin) DHDIR_HANDLE(x,0);
cDhv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
F1HANDLE=@(x,varargin) DDIRH_HANDLE(x,0);
cDvh=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
F1HANDLE=@(x,varargin) HDDIR_HANDLE(x,0);
chDv=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
% 
% F1HANDLE=@(x,varargin) ID2H_HANDLE(x,0);
% cID2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) DI2H_HANDLE(x,0);
cDI2=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

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

ErrC=zeros(NREFS,10);
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
 
   % SP = SparsePattern(NewMesh);
   % SP=1-(SP>0);
    
    M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3());
  %  M_LFE = assemMat_LFE(NewMesh,@MASS_Lump_LFE);

    % Interpolants of test and trial function
    w = assemCochain_1f(NewMesh,F2HANDLE,gauleg(0,1,10));
    %ws = assemCochain_0f(NewMesh,F2sHANDLE);
    u= assemCochain_1f(NewMesh,H_HANDLE,gauleg(0,1,10));
    %Duv= assemCochain_1f(NewMesh,DDIRH_HANDLE,gauleg(0,1,10));

    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation

    % contraction operators
    ContrOne=assemMat_Contr1f(NewMesh,Dir_HANDLE);  % contraction of one forms
    V=Dir_HANDLE(NewMesh.Coordinates);
    %ContrTwo=assemMat_Contr2f_fast(NewMesh,Dir_HANDLE);                 % contraction of two forms
    ContrTwo=assemMat_Contr2f(NewMesh,V);                 % contraction of two forms
    %ContrTwo=assemMat_Contr_1FVertUP(NewMesh,Dir_HANDLE);
    % Contraction Extrusion
    
    % exterior calculus
    ID1 = M*ContrTwo*TopRot;          % (-v x curl u,w)
    %ID1=assemMat_W1F(NewMesh,@STIMA_ContrRot_UP,Dir_HANDLE);
    DI1 = M*TopGrad*ContrOne;         % (grad(v.u),w)
 %   DI2 = -ID1';                                 % (-curl(v x u),w)=(u,v x curl w)
 %   ID2 = -DI1';                                 % (v div u,w)=-(u,grad(v.w))

    LieCompact = assemMat_W1F(NewMesh,@STIMA_LieCompact,Dir_HANDLE);
   
    % Mesh Width
    h(i) = get_MeshWidth(NewMesh);

%    L1c = cID1+cDI1;
    L1c = cDhv+cDvh+chDv;
    L1m = 0.5*cID1+cDvh+chDv;
    
% plot_Norm_W1f(M*LieCompact*u,NewMesh); colorbar;
   % plot_Norm_W1f((ID1+DI1)*u,NewMesh); colorbar;
    % convection discretizations
    [w'*M*LieCompact*u L1c L1m]
    ErrC(i,2) = abs(w'*DI1*u-cDI1);
    ErrC(i,3) = abs(w'*ID1*u-cID1);
    %ErrC(i,3) = abs(w'*DI2*u-cDI2);
    %ErrC(i,4) = abs(w'*ID2*u-cID2);
    %ErrC(i,5) = abs(w'*DI1_FEM*u-cDI1);
    %ErrC(i,6) = abs(w'*ID1_FEM*u-cID1);
    %ErrC(i,7) = abs(w'*DI2_FEM*u-cDI2);
    %ErrC(i,8) = abs(w'*ID2_FEM*u-cID2);
    ErrC(i,9) = abs(w'*(M*LieCompact)*u-L1c);
    ErrC(i,10) = abs(w'*(M*LieCompact)*u-L1m);
    %ErrC(i,10) =abs(w'*(LieInn_FEM+GC_FEM+I        D1_FEM)*u-L1c);
end
% fig = figure('Name','Consistency error convective bi-forms');
% plot(h,h*ErrC(1,2)/h(1),'r-',...
%     h,ErrC(:,2),'ro-',...
%     h,ErrC(:,3),'b-',...
%     h,ErrC(:,4),'bo-',...
%     h,ErrC(:,5),'m-',...
%     h,ErrC(:,6),'mo-',...
%     h,ErrC(:,7),'k-',...
%     h,ErrC(:,8),'ko-',...
%     h,ErrC(:,9),'g-',...
%     h,ErrC(:,10),'go-'); grid('on');
% set(gca,'YScale','log','XScale','log');
% %xlabel('{\bf \Delta t,h= ',num2str(h),'}');
% xlabel(['{\bf h}']);
% ylabel('{\bf L^2-Error}');
% add_Slope(gca,'NorthEast',1);
% add_Slope(gca,'SouthEast',2);
% legend('DI_1','ID_1','DI_2','ID_2','DI_1-FEM','ID_1-FEM','DI_2-FEM','ID_2-FEM','L1c-FEM_p','L1m-FEM_p');

fig = figure('Name','Consistency error convective bi-forms');
plot(h,h*ErrC(1,2)/h(1),'r-',...
    h,ErrC(:,9),'g-',...
    h,ErrC(:,10),'bo-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf Error}');
% add_Slope(gca,'NorthEast',1);
% add_Slope(gca,'SouthEast',2);
legend('h','complete','modified');

