function main_SemiLagrangianSchemesCONSISTENCYLimit(NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
%clear all;b
NREFS_init =0;     % Number of uniform red refinements
% NREFS = 3;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

Dir_Handle=@(x,varargin)straight_flow(x,[0.0,0.0],1);
JvDir_Handle=@(x,varargin)[zeros(size(x)) zeros(size(x))];

H_Handle=@(x,flag,varargin)[ones(size(x,1),1) 0*ones(size(x,1),1)];

F_HANDLE=@(x,flag,varargin) zeros(size(x));

% Load mesh from file
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

err=zeros(4,NREFS);

% refine
for i=1:NREFS_init
    Mesh=refine_REG(Mesh);
end

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
    
    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
    %M = assemMat_Mass1fD(NewMesh);
    
    % discrete pullback
    % Interpolation -> Contraction/Extrusion
    TopGrad=assemMat_TopGrad(NewMesh);   % topological Gradient
    TopRot=assemMat_TopRot(NewMesh);      % topological Rotation
    
    % direct
    ContrOne=assemMat_Contr1f(NewMesh,Dir_Handle);  % contraction of one forms
    ContrTwo=assemMat_Contr2f(NewMesh,Dir_Handle(NewMesh.Coordinates));   % contraction of two form
    Lie_dir = ContrTwo*TopRot+M*TopGrad*ContrOne; ;               % -v x curl u + grad(v.u)
    
    % adjoint
    ContrOne=assemMat_Contr1f(NewMesh,@(x)-Dir_Handle(x));  % contraction of one forms
    ContrTwo=assemMat_Contr2f(NewMesh,-Dir_Handle(NewMesh.Coordinates));   % contraction of two form
    Lie_ad = ContrTwo*TopRot+TopGrad*ContrOne; ;               % -v x curl u + grad(v.u)
    
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
   
    %Dirichlet data
    %system matrix
    A = M;
    %righthand side
    L = M*Lie_dir* H_old;
    %solve system
    LH = A\L;

    %Dirichlet data
    %system matrix
    A = M;
    %righthand side
    %plot_Norm_W1F(Lg,NewMesh); colorbar;
    L = Lie_ad'*M*B_old;
    %solve system
    LB = A\L;

    plot_Norm_W1F(LH,NewMesh);colorbar;
    plot_Norm_W1F(LB,NewMesh);colorbar;
   
    err(1,j)=L2Err_W1F_mod(NewMesh,LH,P7O6(),@(x,varargin)zeros(size(x)),[0.0,0.0],0.5);
    err(2,j)=L2Err_W1F_mod(NewMesh,LB,P7O6(),@(x,varargin)zeros(size(x)),[0.0,0.0],0.5)

    % patches -> projection
    % direct
    Lie_dir = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle,P7O6());
    Lie_dir = Lie_dir+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,Dir_Handle,gauleg(0,1,5));
    Lie_dir = Lie_dir-assemMat_Bnd_LieW1F(NewMesh,[-1],@STIMA_Bnd_LieW1F,Dir_Handle,gauleg(0,1,5));
    Lie_dir = Lie_dir+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,Dir_Handle,gauleg(0,1,5),0.5);
    
    % adjoint
    Lie_ad = assemMat_W1F(NewMesh,@STIMA_ContrRot,@(x,varargin)-Dir_Handle(x),P7O6());
    Lie_ad = Lie_ad+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,5));
    Lie_ad = Lie_ad-assemMat_Bnd_LieW1F(NewMesh,[-1],@STIMA_Bnd_LieW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,5));
    Lie_ad = Lie_ad+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,5),0.5);
        
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    
    %system matrix
    A = M;
    %righthand side
    L = Lie_dir* H_old;
    %solve system
    LH = A\L;

    %system matrix
    A = M;
    %righthand side
    L = Lie_ad'* B_old;
    LB = A\L;

    plot_Norm_W1F(LH,NewMesh);colorbar;
    plot_Norm_W1F(LB,NewMesh);colorbar;
    
    err(3,j)=L2Err_W1F_mod(NewMesh,LH,P7O6(),@(x)zeros(size(x)),[0.0,0.0],0.5);
    err(4,j)=L2Err_W1F_mod(NewMesh,LB,P7O6(),@(x)zeros(size(x)),[0.0,0.0],0.5)
end%  Mesh refinement

figure;
hold on;
plot(mw,err(1,:),'r-x',mw,err(2,:),'g-*',...
    mw,err(3,:),'r--x',mw,err(4,:),'g--*','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('direct Inter','adjoint Inter','direct Patch','adjoint Patch','Location','Northwest');