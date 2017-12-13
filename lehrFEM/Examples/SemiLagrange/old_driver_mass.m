%  driver for Semi-Lagrange, explizit, operor split

%   Copyright 2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

%close all;
clear Mesh;
clear all;
NREFS = 3;
JIG =1;

%
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
%

CFL =[0.8:0.01:1];
CFL =[0.6 0.7 0.8];
T0 = 0;
T1= 0.5;
diff =1;

h = zeros(size(CFL,2),1);
div = zeros(size(CFL,2),3);

v1=1;
v2=1;
Dir_Handle=@(x,t)straight_flow(x);%circ_flow(x);%
F_HANDLE=@(x,flag,t,varargin)zeros(size(x,1),2);
H_Handle=@(x,flag,t,varargin)...
    [pulse_2D(x) pulse_2D(x)];
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
    ContrOne=assemMat_Contr1f(NewMesh,@(x)Dir_Handle(x));  % contraction of one forms
    V=Dir_Handle(NewMesh.Coordinates);
    %ContrTwo=assemMat_Contr2f(NewMesh,V);   % contraction of two forms
    ContrTwo=assemMat_Contr2f_fast(NewMesh,Dir_Handle);   % contraction of two forms

    ID = M*ContrTwo*TopRot;              % -v x curl u geom.
    DI = M*TopGrad*ContrOne;             % grad(v.u) geom.
    
    % pullback of edges
    directions = Dir_Handle(NewMesh.Coordinates,0);
    pbV = trace_vertices(NewMesh,-h(k)*directions);
    
    P = assemMat_SemiLag_W1F(NewMesh, pbV);
        
    % stiffness matrices
    ASL = M + h(k)*diff*C;
    AIE = M + h(k)*diff*C + h(k) * (ID+DI);
    BIE = M ;
    AOS = M + h(k)*diff*C;
    BOS = M - h(k)*(ID+DI);
  
    % timestepping
    H_init = assemLoad_W1F(NewMesh,P7O6(),H_Handle,0);
    HSL_old = M\H_init;
    HIE_old = HSL_old;
    HOS_old = HSL_old;
    
    m_SL = HSL_old'*M*HSL_old;
    m_IE = HIE_old'*M*HIE_old;
    m_OS = HOS_old'*M*HOS_old;
    
    for i = 1:nsteps
        [i i*h(k) nsteps]
        
        [HSL_new,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,1),i*h(k));
        HIE_new = HSL_new;
        HOS_new = HSL_new;

        L = assemLoad_W1F(NewMesh, P7O6(), F_HANDLE,(i)*h(k));
        
        LSL = h(k)*L+M*P*HSL_old-ASL*HSL_new;
        LIE = h(k)*L+BIE*HIE_old-AIE*HIE_new;
        LOS = h(k)*L+BOS*HOS_old-AOS*HOS_new;
        
        HSL_new(FreeDofs) = ASL(FreeDofs,FreeDofs)\LSL(FreeDofs);
        HIE_new(FreeDofs) = AIE(FreeDofs,FreeDofs)\LIE(FreeDofs);
        HOS_new(FreeDofs) = AOS(FreeDofs,FreeDofs)\LOS(FreeDofs);
        
        % Update vectors
        HSL_old=HSL_new;
        HIE_old=HIE_new;
        HOS_old=HOS_new;
        
    end % timestep
    [Dummy,FD_LFE] = assemDir_LFE(NewMesh,-1,@(x,varargin)zeros(size(x,1),1));
    div(k,1) = HSL_old'*M*HSL_old/m_SL;
    div(k,2) = HIE_old'*M*HIE_old/m_IE;
    div(k,3) = HOS_old'*M*HOS_old/m_OS;

end%  Mesh refinement

figure;
plot(h,div(:,1),'bx-', h,div(:,2),'ro-',h,div(:,3),'gd-');
grid('on');
set(gca,'YScale','log');
%set(gca,'YScale','log','XScale','log');
xlabel('\bf \Delta t ');
ylabel('\bf discrete mass');
legend('location','Northwest','Semi-Lagrangian','Explizit Euler','Operator Splitting');
% p = polyfit(log(h),log(Err(:,1)),1);
% add_Slope(gca,'East',1);
