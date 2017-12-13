function main_SemiLagrangianSchemesCONSISTENCY_TIME(CFL,NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
%clear all;b
NREFS_init =1;     % Number of uniform red refinements
% NREFS = 3;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;
v1=0.5;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)v2*(1-x(:,1));
% D1DIR2=@(x)-v2*ones(size(x,1),1);
% D2DIR2=@(x)0*(x(:,2));

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];

% 
% Dir_Handle=@(x,t)straight_flow(x,[0.0,0.0],0.8);
% JvDir_Handle=@(x,t)[zeros(size(x)) zeros(size(x))];

H_Handle=@(x,flag,t,varargin)[ones(size(x,1),1) 0*ones(size(x,1),1)];

F_HANDLE=@(x,flag,t,varargin) zeros(size(x));

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
mw = get_MeshWidth(NewMesh);
NewMesh = init_LEB(NewMesh);
NewMesh=add_Edge2Elem(NewMesh);
NewMesh=add_Patches(NewMesh);
NewMesh = add_DGData(NewMesh);

Lie_dir = assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle,P1O2());
Lie_dir = Lie_dir+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,Dir_Handle,gauleg(0,1,1));
Lie_dir = Lie_dir+assemMat_Bnd_LieW1F(NewMesh,[-1],@STIMA_Bnd_LieW1F,Dir_Handle,gauleg(0,1,5));
Lie_dir = Lie_dir+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,Dir_Handle,gauleg(0,1,1),0.5);

Lie_ad = assemMat_W1F(NewMesh,@STIMA_ContrRot,@(x,varargin)-Dir_Handle(x),P1O2());
Lie_ad = Lie_ad+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LieW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,1));
Lie_ad = Lie_ad+assemMat_Bnd_LieW1F(NewMesh,[-1],@STIMA_Bnd_LieW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,5));
Lie_ad = Lie_ad+assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_LiePennW1F,@(x,varargin)-Dir_Handle(x),gauleg(0,1,1),0.5);

figure; imagesc(Lie_dir); colorbar; title('lim L2 dir');view([0,90])
figure; imagesc(Lie_ad'); colorbar; title('lim L2 ad');view([0,90])

for j=1:NREFS
    h(j) = CFL*mw*2^(-(2+j));
       
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
    %M = assemMat_Mass1fD(NewMesh);
    
    % discrete pullback
    pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h(j));
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h(j));

    % patches 
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    P_p_dir = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    M_q_dir = assemMat_MASSQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    figure; imagesc((M_q_dir-P_p_dir)/h(j)); colorbar; title(' L2 dir');view([0,90])
    err(3,j)=norm(full((M_q_dir-P_p_dir)/h(j)-Lie_dir));
    
    defMesh = NewMesh;
    defMesh.Coordinates = pbVm(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVm);
    P_p_ad = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    M_q_ad = assemMat_MASSQuad_W1F_patches(NewMesh, defMesh,intersec);
    
    figure; imagesc((M_q_ad-P_p_ad')/h(j)); colorbar; title(' L2 ad');view([0,90])
    err(4,j)=norm(full((M_q_ad-P_p_ad')/h(j)-Lie_ad'));
end%  Mesh refinement

figure;
hold on;
plot(h,err(1,:),'r-',h,err(2,:),'g-',...
   h,err(3,:),'r--',h,err(4,:),'g--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('direct Inter','adjoint Inter','direct Patch','adjoint Patch','Location','Northwest');