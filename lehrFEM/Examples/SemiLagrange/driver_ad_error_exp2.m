function driver_ad_error_exp2(CFL)
%  driver for direct Semi-Lagrange (Interpolation and Quadrature) and
%  Eulerian (implicit, semi-implicit, explicit)

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
%clear all;
NREFS =7;
JIG =1;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

%CFL =0.8
T0 = 0;
T1= 0.25;
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
 
% v1=1;
% v2=0.66;
% DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
% D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
% D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);

v1=1;
v2=1;
DIR1=@(x)sin(pi*x(:,1)).*(1-x(:,2).^2);
D1DIR1=@(x)pi*cos(pi*x(:,1)).*(1-x(:,2).^2);
D2DIR1=@(x)sin(pi*x(:,1)).*(-2*x(:,2));
DIR2=@(x)sin(pi*x(:,2)).*(1-x(:,1).^2);
D1DIR2=@(x)sin(pi*x(:,2)).*(-2*x(:,1));
D2DIR2=@(x)pi*cos(pi*x(:,2)).*(1-x(:,1).^2);

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

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
    [ DIR1(x).*(D1H1(x)+D2H2(x)) ...
      DIR2(x).*(D1H1(x)+D2H2(x))];

% -curl(v x h) 
DIH_HANDLE=@(x,flag,t,varargin)T_Handle(t).*...
    [-D2DIR1(x).*H2(x)-DIR1(x).*D2H2(x)+D2DIR2(x).*H1(x)+DIR2(x).*D2H1(x) ...
     +D1DIR1(x).*H2(x)+DIR1(x).*D1H2(x)-D1DIR2(x).*H1(x)-DIR2(x).*D1H1(x)];

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
time1=zeros(NREFS,1);
time2=zeros(NREFS,1);
time3=zeros(NREFS,1);


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

    h = CFL*mw(j)/norm([v1,v2]);
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
    tic;
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());
    toc;
    time3(j)=toc;

    % topological derivatives
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);     % topological Rotation

    % contraction operators
    tic;
    ContrOne=assemMat_Contr1f(NewMesh,@(x)-Dir_Handle(x));  % contraction of one forms
    toc;
    time2(j)=toc;
    %V=Dir_Handle(NewMesh.Coordinates);
    %ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms
    ContrTwo=assemMat_Contr2f_fast(NewMesh,@(x)-Dir_Handle(x));   % contraction of two forms

    ID = -M*ContrTwo*TopRot;              % -v x curl u geom.
    DI = - M*TopGrad*ContrOne;             % grad(v.u) geom.

    % pullback of edges
    % interpolation based
    directions = Dir_Handle(NewMesh.Coordinates,0);
    tic;
    pbVp = trace_vertices(NewMesh,h*directions);
    P = assemMat_SemiLag_W1F(NewMesh, pbVp);
    toc;
    time1(j)=toc;
    % quadrature based
%     defMesh = NewMesh;
%     pbVm = trace_vertices(NewMesh,-h*directions);
%     defMesh.Coordinates = pbVm(:,[1 2]);
%     intersec = aff_elems2(NewMesh, defMesh,pbVm);
%     PV = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
     PV = M*P;

    LieInnM = -assemMat_Lie_Inn_W1F(NewMesh,@(x,t)-Dir_Handle(x,t));
    LieInnP = assemMat_Lie_Inn_W1F(NewMesh,@(x,t)Dir_Handle(x,t));
    LieInn=1/2*(LieInnM+LieInnP);
    LieVol= assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle, P7O6());
    LieVol = LieVol+ assemMat_W1F(NewMesh,@STIMA_GradContr,Dir_Handle, gauleg(0,1,3));
    
   % LieVol = assemMat_W1F(NewMesh,@STIMA_Lie,Dir_Handle,JvDir_Handle);

    % stiffness matrices
    ASL = M + h*diffu*C;
    BSL = P'*M;
    BSLV = PV';
    AIE = M + h*diffu*C - h * (ID' + DI');
    BIE = M;
    AEE = M;
    BEE = M - h*diffu*C + h * (ID' + DI');
    AOS = M + h*diffu*C;
    BOS = M + h*(ID' + DI');
    ALie = M + h*diffu*C - h * (LieInn'+LieVol');
    BLie = M;
    
    % timestepping
    %H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
    H_init = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    HSL_old = H_init;
    HSLV_old =H_init;
    HIE_old = HSL_old;
    HEE_old = HSL_old;
    HOS_old = HSL_old;
    HLie_old = HSL_old;
    time(1)=0;

    for i = 1:nsteps
        [j i i*h nsteps]
        time(i+1)=time(i)+h;

        [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h);
        HSLV_new = HSL_new;
        HIE_new = HSL_new;
        HEE_new = HSL_new;
        HOS_new = HSL_new;
        HLie_new = HSL_new;

        %L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h);
        L = assemCochain_1f(NewMesh, F_HANDLE, gauleg(0,1,10),(i)*h);
        Lee = assemCochain_1f(NewMesh, F_HANDLE, gauleg(0,1,10),(i-1)*h);

        LSL = h*M*L+BSL*HSL_old-ASL*HSL_new;
        LSLV = h*M*L+BSLV*HSLV_old-ASL*HSLV_new;
        LIE = h*M*L+BIE*HIE_old-AIE*HIE_new;
        LEE = h*M*Lee+BEE*HEE_old-AEE*HEE_new;
        LOS = h*M*L+BOS*HOS_old-AOS*HOS_new;
        LLie = h*M*L+BLie*HLie_old-ALie*HLie_new;
        tic
        HSL_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSL(FreeDofs);
        toc; time1(j)=time1(j)+toc;
        %HSLV_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSLV(FreeDofs);
        tic;
        HIE_new(FreeDofs) = AIE(FreeDofs,FreeDofs)\LIE(FreeDofs);
        toc; time1(j)=time1(j)+toc;
        %HEE_new(FreeDofs) = AEE(FreeDofs,FreeDofs)\LEE(FreeDofs);
        HOS_new(FreeDofs) = AOS(FreeDofs,FreeDofs)\LOS(FreeDofs);
        HLie_new(FreeDofs) = ALie(FreeDofs,FreeDofs)\LLie(FreeDofs);

        % Update vectors
        HSL_old=HSL_new;
        HSLV_old=HSLV_new;
        HIE_old=HIE_new;
        HEE_old=HEE_new;
        HOS_old=HOS_new;
        HLie_old=HLie_new;
    end % timestep
    
    Err(j,1) = L2Err_W1F(NewMesh,HSL_old,P7O6(),H_Handle,0,i*h);
    Err(j,2) = L2Err_W1F(NewMesh,HIE_old,P7O6(),H_Handle,0,i*h);
    Err(j,3) = L2Err_W1F(NewMesh,HOS_old,P7O6(),H_Handle,0,i*h);
    Err(j,4) = L2Err_W1F(NewMesh,HEE_old,P7O6(),H_Handle,0,i*h);
    Err(j,5) = L2Err_W1F(NewMesh,HSLV_old,P7O6(),H_Handle,0,i*h);
    Err(j,6) = L2Err_W1F(NewMesh,HLie_old,P7O6(),H_Handle,0,i*h);

end%  Mesh refinement


filename=['./figures/L2ErrAdjoint2',sprintf('%1.2f',CFL),'.mat'];
save(filename, 'Err','mw');

figure();
%plot(CFL,Err(:,1),'bx-', CFL,Err(:,2),'ro-',CFL,Err(:,3),'gd-',CFL,Err(:,6),'c>-',CFL,Err(:,5),'y<-',CFL,Err(:,4),'kx-');
plot(mw,mw*Err(1,1)/mw(1),'k--',mw,Err(:,1),'bx-', mw,Err(:,2),'ro-',mw,Err(:,3),'gd-',mw,Err(:,6),'c>-','linewidth', 3);
grid('on');
title('Adjoint method');
set(gca,'YScale','log','XScale','log','FontWeight','bold');
%set(gca,'YScale','log','XScale','log');
xlabel('\bf h');
ylabel('\bf L^2-error');
%legend('location','Northwest','Semi-Lagrangian','Eulerian, implicit','Eulerian, semi-implicit','Eulerian, explicit','Semi-Lagrangian_{var}');
legend('location','Northwest','h','\bf Semi-Lagrangian','\bf Eulerian Upwind, implicit','\bf Eulerian Upwind, semi-implicit','\bf Eulerian Standard, implicit');
% p = polyfit(log(h),log(Err(:,1)),1);
% add_Slope(gca,'East',1);
filename=['./figures/L2ErrAdjoint2',sprintf('%1.2f',CFL),'.fig'];
saveas(gcf,filename);
filename=['./figures/L2ErrAdjoint2',sprintf('%1.2f',CFL),'.eps'];
print('-depsc', filename);
