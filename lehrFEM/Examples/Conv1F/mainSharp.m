% Run script discrete differential forms.

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

clear Mesh;
tic;
% close all;
% clear all;

onlyStandard = 1;

% Initialize constant
JIG=1;
NREFS =5;
Sigma = 0;
sharp = 0;

EPSI_Handle = @(x,varargin)ones(size(x,1),1);
MU_HANDLE=@(x,varargin)1;

H1=@(x) x(:,1)-x(:,1).^2;
H2=@(x) x(:,1)-x(:,1).^2;

DIR1 = @(x)  0*ones(size(x(:,1)));
D1DIR1 = @(x)  0*ones(size(x(:,1)));
D2DIR1 = @(x)  0*ones(size(x(:,1)));
DIR2=@(x) 1*ones(size(x(:,1)));
D1DIR2=@(x) 0*ones(size(x(:,1)));
D2DIR2=@(x) 0*ones(size(x(:,1)));

V_Handle=@(x,varargin)[DIR1(x), DIR2(x)];
JacV_Handle=@(x,varargin)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
U_EX_Handle=@(x,varargin)[H1(x) H2(x)];

QuadRule2D=gauleg(0,1,10);
QuadRule=Duffy(TProd(gauleg(0,1,10)));

L2err=zeros(NREFS,1);
h=zeros(NREFS,1);
Dofs=zeros(NREFS,1);
approx=zeros(NREFS,3);

if ~sharp
    % Load mesh from file
    Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

    % Add edge data structure
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    % Mesh.BdFlags(Loc) = [-2,-2,-2,-1,-2,-1,-1,-1];
    Mesh.BdFlags(Loc) = [-1,-2,-1,-2,-2,-2,-2,-2];
end

for i = 1:NREFS

    if sharp
        Mesh=create_Mesh(i,Sigma);
        Mesh= setBndFlags(Mesh);
        Mesh = orient_Elems(Mesh);
    else
         Mesh=refine_REG(Mesh);
    end
        
    % Mesh preprocessing
    switch(JIG)
        case 1
            New_Mesh = Mesh;
        case 2
            Loc = get_BdEdges(Mesh);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            FixedPos = zeros(size(Mesh.Coordinates,1),1);
            FixedPos(Loc) = 1;
            New_Mesh = jiggle(Mesh,FixedPos);
        case 3
            Loc = get_BdEdges(Mesh);
            Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
            FixedPos = zeros(size(Mesh.Coordinates,1),1);
            FixedPos(Loc) = 1;
            New_Mesh = smooth(Mesh,FixedPos);
    end
    
    New_Mesh=add_Edge2Elem(New_Mesh); 
    New_Mesh = add_DGData(New_Mesh);

    % stiffness matrices
    M_W1F = assemMat_W1F(New_Mesh,@MASS_W1F,MU_HANDLE, P7O6());  % Mass Matrix Whitney-1-Forms
    
    % boundary data on inflow boundary and righthandside
    L=assemCochain_1f(New_Mesh,@(x)zeros(size(x)),gauleg(0,1,10));     % -v x curl u righthand side
    
    Lie1 = assemMat_W1F(New_Mesh,@STIMA_ContrRot,V_Handle,P7O6());
    Lie2 =  assemMat_Inn_LieW1F(New_Mesh,@STIMA_Inn_LieW1F,V_Handle,QuadRule2D);
    Lie2bnd = assemMat_Bnd_LieW1F(New_Mesh,[-2],@STIMA_Bnd_LieW1F,V_Handle,gauleg(0,1,5));
    Lie2pen = assemMat_Inn_LieW1F(New_Mesh,@STIMA_Inn_LiePennW1F,V_Handle,QuadRule2D,0.5);
    
    % exact solution
    cU=assemCochain_1f(New_Mesh,U_EX_Handle,gauleg(0,1,10));
        
    % calculate solutions
    % free DOFS
    nonInFlow=find(New_Mesh.BdFlags ~= -1);
    nonOutFlow=find(New_Mesh.BdFlags ~= -2);
    FreeDofs=nonInFlow;
    
    % penalized standard scheme
    S_pst=Lie1+Lie2+Lie2pen-Lie2bnd;
    cfbnd = zeros(size(L));
    cfbnd = assemLoad_Bnd_LieW1F(New_Mesh,[-1],cfbnd,gauleg(0,1,10),U_EX_Handle, V_Handle);
    cU(FreeDofs) = 0;
    Lf_pst = L+cfbnd-S_pst*cU;

    % solve system
    cU(FreeDofs) = S_pst(FreeDofs,FreeDofs)\Lf_pst(FreeDofs);
    
    %plot_Norm_W1F(cU,New_Mesh); colorbar; title('U Patches');
    
    L2err(i) = L2Err_W1F(New_Mesh,cU,QuadRule,U_EX_Handle)
    
    h(i) = get_MeshWidth(New_Mesh);
    Dofs(i) = size(New_Mesh.Edges,1);
end

slope = [diff(log(L2err))./diff(log(h))]

fig = figure('Name','L^2-error');
plot(h,L2err(1)/h(1)*h,'k--',h,L2err(1)/h(1)^(1/2)*h.^(1/2),'k--',h,L2err,'r-'); grid('on');
set(gca,'XScale','log','YScale','log');
xlabel('{\bf h}');
ylabel('{\bf L^2-Error}');
legend('h','h^{1/2}','stabelized','Location','SouthEast')