% consistency of bilinarforms approximating the LieDerivative 
% parts i_v d und d i_v for 1 and 2 forms, using extrusion contraction   

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all
NREFS =7;  % Number of unifrom red refinements
JIG =2;

H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);

H_HANDLE=@(x,flag,varargin)...
    [ H1(x) H2(x)];

B_HANDLE=@(x,varargin)[D1H2(x)-D2H1(x)];

v1=-0.66;
v2=-0.3;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_HANDLE=@(x)[DIR1(x), DIR2(x)];

% -(v x b) 
IB_HANDLE=@(x,flag,varargin)[-DIR2(x).*B_HANDLE(x) DIR1(x).*B_HANDLE(x)];

% -v x curl h
ID1H_HANDLE=@(x,flag,varargin)...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
       DIR1(x).*(D1H2(x)-D2H1(x))];
  
% exact values of bi-forms 
F2HANDLE=@(x,varargin) [x(:,1).^2+x(:,2) x(:,2).*x(:,1)];

F1HANDLE=@(x,varargin) IB_HANDLE(x,0);
cIB=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) ID1H_HANDLE(x,0);
cID1=innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

% Load mesh from file

% initial mesh
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');

% Add edge data structure

Mesh = add_Edges(Mesh);
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

Err=zeros(NREFS,4);
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
    M_P = assemMat_P0(NewMesh,@MASS_P0);
    
    % Interpolants of test and trial function
    W = assemCochain_1F(NewMesh,F2HANDLE,gauleg(0,1,10));
    w = W;
    H = assemCochain_1F(NewMesh,H_HANDLE,gauleg(0,1,10));
    H = H;

    u = assemCochain_2F(NewMesh,B_HANDLE,P7O6);
    %u=M_P\U;
    
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation
    
    % contraction operators
    V=Dir_HANDLE(NewMesh.Coordinates);
    ContrTwo2=assemMat_Contr2f_fast(NewMesh,Dir_HANDLE);                 % contraction of two forms
    ContrTwo=assemMat_Contr2f(NewMesh,V);                 % contraction of two forms
 
    % Mesh Width
    h(i) = get_MeshWidth(NewMesh);
    
    Err(i,1) = abs(w'*M*ContrTwo*u-cIB);
    Err(i,2) = abs(w'*M*ContrTwo*TopRot*H-cID1);
    Err(i,3) = abs(w'*M*ContrTwo2*TopRot*H-cID1)
    
end
fig = figure('Name','Consistency error bi-forms');
plot(h,Err(:,1),'r+-',h,Err(:,2),'b+-',h,Err(:,3),'g+-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf Error}');
add_Slope(gca,'NorthEast',1);
legend('IB','IDH','Interpol');

% [diff(log(h)).\diff(log(Err(:,1)))...
%     diff(log(h)).\diff(log(Err(:,2)))...
%     diff(log(h)).\diff(log(Err(:,3)))...
%     diff(log(h)).\diff(log(Err(:,4)))...
%     diff(log(h)).\diff(log(Err(:,5)))...
%     diff(log(h)).\diff(log(Err(:,6)))...
%     diff(log(h)).\diff(log(Err(:,7)))...
%     diff(log(h)).\diff(log(Err(:,8)))]