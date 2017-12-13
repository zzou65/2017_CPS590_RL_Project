%  driver for Semi-Lagrange, implizit, operor split

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS =3;
JIG =2;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

%CFL =[0.8:0.01:1];
%phi=1;
CFL =[0.6 0.7 0.8 0.9 1];
%CFL =[0.1 0.2 0.3 0.4];
T0 = 0;
T1= 2;
diff =1 ;

h = zeros(size(CFL,2),1);
Err = zeros(size(CFL,2),4);
div = zeros(size(CFL,2),4);

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
% 
% % div free
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
% 
% H1 = @(x) 4 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
% D1H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
% D2H1=@(x) 4 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
% D11H1=@(x) -16 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
% D12H1=@(x) -16 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
% D22H1=@(x) -24 * (1-x(:,1).^2).^2 .*x(:,2);
% H2=@(x) -4 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
% D1H2=@(x) -4 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
% D2H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
% D11H2=@(x) 24 * x(:,1).*(1-x(:,2).^2).^2;
% D12H2=@(x) 16 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
% D22H2=@(x) 16 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);

v1=0.66;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

% H1=@(x)sin(pi.*x(:,1)).*(1-x(:,2));
% D1H1=@(x)pi*(1-x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)-sin(pi.*x(:,1));
% D11H1=@(x)-pi^2*(1-x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)-pi.*cos(pi.*x(:,1));
% D22H1=@(x)0*x(:,1);
% 
% H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
% D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
% D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
% D11H2=@(x)(1-x(:,2).^2).*(-2);
% D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
% D22H2=@(x)(-2).*(1-x(:,1).^2);
% 
% v1=0.66;
% v2=1;
% DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

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
    diff * CCH_HANDLE(x,flag,t) + ...
    DtH_HANDLE(x,flag,t) + ...
    IDH_HANDLE(x,flag,t) + ...
    DIH_HANDLE(x,flag,t);

G_HANDLE=@(x,flag,t,varargin)...
    diff * CCH_HANDLE(x,flag,t) + ...
    DtH_HANDLE(x,flag,t) + ...
    IDH_HANDLE(x,flag,t) + ...
    DIH_HANDLE(x,flag,t);


% DIR1=@(x)0.25*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% DIR2=@(x)1*sin(pi.*x(:,2)).*sin(pi.*x(:,1));%(1-x(:,2).^2).*(1-x(:,1).^2);
Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),2);

% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

% refine
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
    
end

for k = 1:size(CFL,2)
    
    mw = get_MeshWidth(NewMesh);

    h(k) = CFL(k)*mw/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h(k));
    h(k) = h(k)-(T0+nsteps*h(k)-T1)/nsteps;
    steps = nsteps;

    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);

    % plot_Mesh(NewMesh,'petas');

    nEdges= size(NewMesh.Edges,1);

    % Assemble Curl-curl matrix, MASS matrix and load vector
    C = assemMat_W1F(NewMesh,@STIMA_Curl_W1F,SIGMA_HANDLE,P7O6());
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());
    
    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);     % topological Rotation
     
    % contraction operators
    ContrOne=assemMat_Contr1f(NewMesh,@(x)-Dir_Handle(x));  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    ContrTwo=assemMat_Contr2f(NewMesh,-V);   % contraction of two forms

    ID = -M*ContrTwo*TopRot;              % -v x curl u geom.
    DI = -M*TopGrad*ContrOne;             % grad(v.u) geom.
    
    % pullback of edges
    directions = Dir_Handle(NewMesh.Coordinates,0);
    pbV = trace_vertices(NewMesh,h(k)*directions);
    
    P = assemMat_SemiLag_W1F(NewMesh, pbV);
        
    % stiffness matrices
    ASL = M + h(k)*diff*C;
    AIE = M + h(k)*diff*C - h(k) * (ID'+DI');
    BIE = M ;
    AEE = M;
    BEE = M - h(k)*diff*C + h(k) * (ID'+DI');;
    AOS = M + h(k)*diff*C;
    BOS = M + h(k)*(ID'+DI');
  
    % timestepping
    H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
    HSL_old = M\H_init;
    HIE_old = HSL_old;
    HEE_old = HSL_old;
    HOS_old = HSL_old;
    
    for i = 1:nsteps
        [i i*h(k) nsteps]
        
        [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h(k));
        HIE_new = HSL_new;
        HEE_new = HSL_new;
        HOS_new = HSL_new;

        L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h(k));
        
        LSL = h(k)*L+P'*M*HSL_old-ASL*HSL_new;
        LIE = h(k)*L+BIE*HIE_old-AIE*HIE_new;
        LEE = h(k)*L+BEE*HEE_old-AEE*HEE_new;
        LOS = h(k)*L+BOS*HOS_old-AOS*HOS_new;
        
        HSL_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSL(FreeDofs);
        HIE_new(FreeDofs) = AIE(FreeDofs,FreeDofs)\LIE(FreeDofs);
        HEE_new(FreeDofs) = AEE(FreeDofs,FreeDofs)\LEE(FreeDofs);
        HOS_new(FreeDofs) = AOS(FreeDofs,FreeDofs)\LOS(FreeDofs);
        
        % Update vectors
        HSL_old=HSL_new;
        HIE_old=HIE_new;
        HEE_old=HEE_new;
        HOS_old=HOS_new;
        
    end % timestep
    Err(k,1) = L2Err_W1F(NewMesh,HSL_old,P7O6(),H_Handle,0,(i)*h(k));
    Err(k,2) = L2Err_W1F(NewMesh,HIE_old,P7O6(),H_Handle,0,(i)*h(k));
    Err(k,3) = L2Err_W1F(NewMesh,HOS_old,P7O6(),H_Handle,0,(i)*h(k));
    Err(k,4) = L2Err_W1F(NewMesh,HEE_old,P7O6(),H_Handle,0,(i)*h(k));
    
    [Dummy,FD_LFE] = assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
    div(k,1) = norm(TopGrad(:,FD_LFE)'*M*HSL_old);
    div(k,2) = norm(TopGrad(:,FD_LFE)'*M*HIE_old);
    div(k,3) = norm(TopGrad(:,FD_LFE)'*M*HOS_old);
    div(k,4) = norm(TopGrad(:,FD_LFE)'*M*HEE_old);

end%  Mesh refinement

figure('Name','Error');
plot(CFL,Err(:,1),'bx-', CFL,Err(:,2),'ro-',CFL,Err(:,3),'gd-',CFL,Err(:,4),'cd-');
grid('on');
set(gca,'YScale','log');
%set(gca,'YScale','log','XScale','log');
xlabel('CFL ');
legend('location','Northwest','Semi-Lagrangian','Implizit Euler','Operator Splitting','Explizit Euler');
% p = polyfit(log(h),log(Err(:,1)),1);
% add_Slope(gca,'East',1);

figure('Name','weak divergence')
plot(h,div(:,1),'bx-', h,div(:,2),'ro-',h,div(:,3),'gd-',h,div(:,4),'cd-');
grid('on');
set(gca,'YScale','log');
%set(gca,'YScale','log','XScale','log');
xlabel('\bf \Delta t ');
ylabel('\bf weak divergence');
legend('location','Northwest','Semi-Lagrangian','Implizit Euler','Operator Splitting','Explicit Euler');
% p = polyfit(log(h),log(Err(:,1)),1);
% add_Slope(gca,'East',1);
