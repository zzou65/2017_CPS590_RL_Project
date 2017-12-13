% consistency of bilinarforms approximating the LieDerivative 
% parts i_v d und d i_v for 1 and 2 forms, using extrusion contraction   

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all
NREFS = 5;  % Number of unifrom red refinements
JIG =1;

% H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

H1=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H1=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H1=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% 
H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% 
% H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% % 
%H2=@(x)2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
%D1H2=@(x)2*pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
%D2H2=@(x)2*pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

% H1 = @(x) 4 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2H1=@(x) 4 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% H2=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1H2=@(x) -4 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);

% H1 = @(x) (1-x(:,1).^2);
% D1H1 = @(x) (-2*x(:,1));
% D2H1 = @(x) 0*x(:,1);
% H2 = @(x) (1+2*x(:,2).*x(:,1));
% D1H2 = @(x) (2*x(:,2));
% D2H2 = @(x) (2*x(:,1));

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
%DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
%D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
%D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

Dir_HANDLE=@(x)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

H_HANDLE=@(x,flag,varargin)...
    [H1(x) H2(x)];

%  -v x curl h
ID1H_HANDLE=@(x,flag,varargin)...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
       DIR1(x).*(D1H2(x)-D2H1(x))];
   
%  h x curl v
IDAH_HANDLE=@(x,flag,varargin)...
    [H2(x).*(D1DIR2(x)-D2DIR1(x)) ...
     -H1(x).*(D1DIR2(x)-D2DIR1(x))];

% grad v h
DI1H_HANDLE=@(x,flag,varargin)...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
     D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];
 
 % v * div h
ID2H_HANDLE=@(x,flag,varargin)...
    [DIR1(x).*(D1H1(x)+D2H2(x)) ...
     DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h) 
DI2H_HANDLE=@(x,flag,varargin)...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
     +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];
 
%exact values of bi-forms 
G1=@(x) x(:,1).^2+x(:,2); G2=@(x) x(:,2).*x(:,1);
%G1=@(x) 1-x(:,2); G2=@(x) -1+x(:,1);
%G2=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
%G1=@(x) (1-x(:,1).^2); G2=@(x) (1+2*x(:,2).*x(:,1));

%  g x curl v
DIAH_HANDLE=@(x,flag,varargin)...
    [G2(x).*(D1DIR2(x)-D2DIR1(x)) ...
     -G1(x).*(D1DIR2(x)-D2DIR1(x))];


F2HANDLE=@(x,varargin) [G1(x) G2(x)];

F1HANDLE=@(x,varargin) ID1H_HANDLE(x,0);
cID1 = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE=@(x,varargin) IDAH_HANDLE(x,0);
cIDA = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE = @(x,varargin) DI1H_HANDLE(x,0);
cDI1 = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE = @(x,varargin) ID2H_HANDLE(x,0);
cID2 = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE = @(x,varargin) DI2H_HANDLE(x,0);
cDI2 = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

F1HANDLE = @(x,varargin) DIAH_HANDLE(x,0);
cDIA = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);

[cIDA cID1 cDI1 cID2 cDI2]
% Load mesh from file

% initial mesh
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');

% Add edge data structure

Mesh = add_Edges(Mesh);
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

Err = zeros(NREFS,8);
h = zeros(NREFS,1);

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
    NewMesh = add_Edge2Elem(NewMesh);
    NewMesh = add_Patches(NewMesh);

    M = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P3O3());
    Mb = assemMat_W1F(NewMesh,@MASS_W1F,@(x,varargin)1, P1O2());
    
    % Interpolants of test and trial function
    w = assemCochain_1f(NewMesh,F2HANDLE,gauleg(0,1,10));
    u = assemCochain_1f(NewMesh,H_HANDLE,gauleg(0,1,10));
    V = Dir_HANDLE(NewMesh.Coordinates);
    h(i) = get_MeshWidth(NewMesh);
    
    % Semi-Lagrange
    dt=0.1*h(i)/norm([v1,v2]);
    
    %pb_c = assemCochain_1f(NewMesh,F1HANDLE,gauleg(0,1,10));
    %pb_ad_c = assemCochain_1f(NewMesh,F2HANDLE,gauleg(0,1,10));
   % plot_Norm_W1F(B1b*u,NewMesh); colorbar;
   % plot_Norm_W1F(M*A1*u,NewMesh); colorbar;
 
   % direct
    pulled_back_vertices = trace_vertices(NewMesh,-dt*V);
    A1 = assemMat_SemiLag_W1F(NewMesh, pulled_back_vertices);
    pbB = trace_bcenters(NewMesh, Dir_HANDLE,JvDir_Handle,-dt);
    B1b = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);
    % B1s = assemMat_SemiLagQuad_W1F_strang1(NewMesh,
    % @(x)Dir_HANDLE(x),JvDir_Handle,-dt);
    
    F1HANDLE=@(x,varargin) [(1-dt*D1DIR1(x)).*H1(x-dt*Dir_HANDLE(x))+(-dt*D1DIR2(x)).*H2(x-dt*Dir_HANDLE(x)) ...
                                         (-dt*D2DIR1(x)).*H1(x-dt*Dir_HANDLE(x))+(1-dt*D2DIR2(x)).*H2(x-dt*Dir_HANDLE(x))];
    pb = innerproduct(F1HANDLE,F2HANDLE,-1,1,-1,1);
    Err(i,1) = abs(w'*B1b*u-pb);
    % Err(i,2) = abs(w'*B1s*u-pb); 
    Err(i,2) = L2Err_W1F(NewMesh,M\B1b*u,P7O6(),F1HANDLE,0);
    Err(i,3) = abs(w'*M*A1*u-pb);
    Err(i,4) = L2Err_W1F(NewMesh,A1*u,P7O6(),F1HANDLE,0);
    % Err(i,4) = abs(w'*B1*u-pb);
    
    % adjoint
    pulled_back_vertices = trace_vertices(NewMesh,dt*V);
    A2 = assemMat_SemiLag_W1F(NewMesh, pulled_back_vertices);
    % B2 = assemMat_SemiLagQuad_W1F(NewMesh, pulled_back_vertices,JDir_Handle,-dt);
    pbB = trace_bcenters(NewMesh, Dir_HANDLE,JvDir_Handle,dt);
    B2b = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);
    % B2s = assemMat_SemiLagQuad_W1F_strang1(NewMesh, @(x)Dir_HANDLE(x),JvDir_Handle,dt);
   
    F1HANDLE=@(x,varargin) [(1+dt*D1DIR1(x)).*G1(x+dt*Dir_HANDLE(x))+(dt*D1DIR2(x)).*G2(x+dt*Dir_HANDLE(x)) ...
                                         (dt*D2DIR1(x)).*G1(x+dt*Dir_HANDLE(x))+(1+dt*D2DIR2(x)).*G2(x+dt*Dir_HANDLE(x))];
    pb_ad = innerproduct(@(x)H_HANDLE(x,0),F1HANDLE,-1,1,-1,1);

    Err(i,5) = abs(w'*B2b'*u-pb_ad);
   % Err(i,6) = abs(w'*B2s'*u-pb_ad);
    Err(i,6) = L2Err_W1F(NewMesh,M\B2b*w,P7O6(),F1HANDLE,0);
    Err(i,7) = abs(w'*A2'*M*u-pb_ad);
    Err(i,8) = L2Err_W1F(NewMesh,A2*w,P7O6(),F1HANDLE,0);
   % Err(i,8) = abs(w'*B2'*u-pb_ad);
   
    Err(i,9) = abs((w'*M*u-pb)/dt-(cID1+cDI1));
    Err(i,10) = abs((w'*M*u-pb_ad)/dt-(cID2+cDI2))
    
end
fig = figure('Name','Consistency error bi-forms 1-forms');
plot(h,Err(:,1),'r+-',...
    h,Err(:,2),'go-',...
    h,Err(:,3),'b*-',...
    h,Err(:,4),'c.-',...
    h,Err(:,9),'md-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf Error}');
add_Slope(gca,'NorthEast',1);
legend('bary','strang','Interpol','vertex','timedis');

fig = figure('Name','Consistency error bi-forms 2-forms');
plot(h,Err(:,5),'r+-',...
    h,Err(:,6),'go-',...
    h,Err(:,7),'b*-',...
    h,Err(:,8),'c.-',...
    h,Err(:,10),'md-'); grid('on');
set(gca,'YScale','log','XScale','log');
%xlabel('{\bf \Delta t,h= ',num2str(h),'}');
xlabel(['{\bf h}']);
ylabel('{\bf Error}');
add_Slope(gca,'NorthEast',1);
legend('bary','strang','Interpol','vertex','timedis');

[diff(log(h)).\diff(log(Err(:,1)))...
    diff(log(h)).\diff(log(Err(:,2)))...
    diff(log(h)).\diff(log(Err(:,3)))...
    diff(log(h)).\diff(log(Err(:,4)))...
    diff(log(h)).\diff(log(Err(:,5)))...
    diff(log(h)).\diff(log(Err(:,6)))...
    diff(log(h)).\diff(log(Err(:,7)))...
    diff(log(h)).\diff(log(Err(:,8)))...
    diff(log(h)).\diff(log(Err(:,9)))...
    diff(log(h)).\diff(log(Err(:,10)))]