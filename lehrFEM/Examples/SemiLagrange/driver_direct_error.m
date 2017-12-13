%  driver for direct Semi-Lagrange (Interpolation and Quadrature) and
%  Eulerian (implicit, semi-implicit, explicit)

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS =4;
JIG =1;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL =0.1;
T0 = 0;
T1= 0.5;
diffu =1;

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

% H1=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H1=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H1=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% D11H1=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D12H1=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
% D22H1=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% H2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D1H2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
% D2H2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
% D11H2=@(x)-pi.^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
% D12H2=@(x)pi^2*cos(pi.*x(:,2)).*cos(pi.*x(:,1));
% D22H2=@(x)-pi^2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));

v1=0.5;
v2=0.5;
DIR1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
D1DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
D2DIR1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
D11DIR1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
D12DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
D22DIR1=@(x) -1 * (1-x(:,1).^2).^2 .*x(:,2);
DIR2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
D1DIR2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
D2DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
D11DIR2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
D12DIR2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
D22DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
% 
v1=1;
v2=0.66;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

DIR1=@(x)sin(pi*x(:,1)).*(1-x(:,2).^2);
D1DIR1=@(x)pi*cos(pi*x(:,1)).*(1-x(:,2).^2);
D2DIR1=@(x)sin(pi*x(:,1)).*(-2*x(:,2));
DIR2=@(x)sin(pi*x(:,2)).*(1-x(:,1).^2);
D1DIR2=@(x)sin(pi*x(:,2)).*(-2*x(:,1));
D2DIR2=@(x)pi*cos(pi*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];

H_Handle=@(x,flag,t,varargin)...
    T_Handle(t).*[ H1(x) H2(x)];

% curl curl h
CCH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D22H1(x) + D12H2(x) ...
    D12H1(x)-D11H2(x)];

% d_t h
DtH_HANDLE=@(x,flag,t,varargin)...
    DT_Handle(t).*[ H1(x) H2(x)];

%  -v x curl h
IDH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
    DIR1(x).*(D1H2(x)-D2H1(x))];

% grad v h
DIH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
    D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];

F_HANDLE=@(x,flag,t,varargin)...
    diffu * CCH_HANDLE(x,flag,t) + ...
    DtH_HANDLE(x,flag,t) + ...
    IDH_HANDLE(x,flag,t) + ...
    DIH_HANDLE(x,flag,t);

% Load mesh from file

Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure

Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

div=zeros(NREFS,3);
mw=zeros(NREFS,1);

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

    mw(j) = get_MeshWidth(NewMesh);

    h = CFL*sqrt(mw(j))/norm([v1,v2]);
    nsteps = ceil((T1-T0)/h);
    h = h-(T0+nsteps*h-T1)/nsteps;
    steps = nsteps;

    NewMesh = init_LEB(NewMesh);
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
    ContrOne=assemMat_Contr1f(NewMesh,@(x)Dir_Handle(x));  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    %ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms
    ContrTwo=assemMat_Contr2f_fast(NewMesh,Dir_Handle);   % contraction of two forms

    ID = M*ContrTwo*TopRot;              % -v x curl u geom.
    DI = M*TopGrad*ContrOne;             % grad(v.u) geom.

    % pullback of edges
    % interpolation based
    directions = Dir_Handle(NewMesh.Coordinates,0);
    pbVm = trace_vertices(NewMesh,-h*directions);
    P = assemMat_SemiLag_W1F(NewMesh, pbVm);
    % quadrature based
    %defMesh = NewMesh;
    %pbVp = trace_vertices(NewMesh,h*directions);
    %defMesh.Coordinates = pbVp(:,[1 2]);
    %intersec = aff_elems2(NewMesh, defMesh,pbVp);
    %PV = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    %PV = assemMat_SemiLagQuad_W1F_patches2nd(NewMesh, defMesh,intersec);
    PV = M*P;
    
    % stiffness matrices
    ASL = M + h*diffu*C;
    BSL = M*P;
    BSLV = PV;
    AIE = M + h*diffu*C + h * (ID+DI);
    BIE = M;
    AEE = M;
    BEE = M - h*diffu*C - h * (ID+DI);
    AOS = M + h*diffu*C;
    BOS = M - h*(ID+DI);

    % timestepping
    %H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
    H_init = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    HSL_old = H_init;
    HSLV_old =H_init;
    HIE_old = HSL_old;
    HEE_old = HSL_old;
    HOS_old = HSL_old;

    time(1)=0;

    for i = 1:nsteps
        [i i*h nsteps]
        time(i+1)=time(i)+h;

        [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h);
        HSLV_new = HSL_new;
        HIE_new = HSL_new;
        HEE_new = HSL_new;
        HOS_new = HSL_new;

        %L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h);
        L = assemCochain_1f(NewMesh, F_HANDLE, gauleg(0,1,10),(i)*h);
        Lee = assemCochain_1f(NewMesh, F_HANDLE, gauleg(0,1,10),(i)*h);

        LSL = h*M*L+BSL*HSL_old-ASL*HSL_new;
        LSLV = h*M*L+BSLV*HSLV_old-ASL*HSLV_new;
        LIE = h*M*L+BIE*HIE_old-AIE*HIE_new;
        LEE = h*M*Lee+BEE*HEE_old-AEE*HEE_new;
        LOS = h*M*L+BOS*HOS_old-AOS*HOS_new;

        HSL_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSL(FreeDofs);
        HSLV_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSLV(FreeDofs);
        HIE_new(FreeDofs) = AIE(FreeDofs,FreeDofs)\LIE(FreeDofs);
        HEE_new(FreeDofs) = AEE(FreeDofs,FreeDofs)\LEE(FreeDofs);
        HOS_new(FreeDofs) = AOS(FreeDofs,FreeDofs)\LOS(FreeDofs);

        % Update vectors
        HSL_old=HSL_new;
        HSLV_old=HSLV_new;
        HIE_old=HIE_new;
        HEE_old=HEE_new;
        HOS_old=HOS_new;

    end % timestep
    
    Err(j,1) = L2Err_W1F(NewMesh,HSL_old,P7O6(),H_Handle,0,i*h);
    Err(j,2) = L2Err_W1F(NewMesh,HIE_old,P7O6(),H_Handle,0,i*h);
    Err(j,3) = L2Err_W1F(NewMesh,HOS_old,P7O6(),H_Handle,0,i*h);
    Err(j,4) = L2Err_W1F(NewMesh,HEE_old,P7O6(),H_Handle,0,i*h);
    Err(j,5) = L2Err_W1F(NewMesh,HSLV_old,P7O6(),H_Handle,0,i*h);

end%  Mesh refinement

figure('Name','Error');
plot(mw,mw*(Err(1,1)/mw(1)),'k--',mw,Err(:,1),'bx-', mw,Err(:,2),'ro-',mw,Err(:,3),'gd-',mw,Err(:,4),'c>-',mw,Err(:,5),'y<-');
grid('on');
title('Direct Method');
set(gca,'YScale','log','XScale','log');
%set(gca,'YScale','log','XScale','log');
xlabel('h ');
legend('location','Northwest','h','Semi-Lagrangian','Eulerian, implicit','Eulerian, semi-implicit','Eulerian, explicit','Semi-Lagrangian_{var}');
% p = polyfit(log(h),log(Err(:,1)),1);
% add_Slope(gca,'East',1);
filename=['./figures/L2ErrDirect',sprintf('%1.2f',CFL),'.fig'];
saveas(gcf,filename);
filename=['./figures/L2ErrDirect',sprintf('%1.2f',CFL),'.eps'];
print('-depsc', filename);
filename=['./figures/L2ErrDirect',sprintf('%1.2f',CFL),'.mat'];
save(filename, 'Err','mw');
