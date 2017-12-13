function main_SemiLagrangianSchemesCONSISTENCY(CFL,NREFS)
%  Semi-Lagrange- version for MHD

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%close all;

clear Mesh;
%clear all;b
NREFS_init =0;     % Number of uniform red refinements
TimeNREFs = 4;
% NREFS = 3;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

Dir_Handle=@(x,t)straight_flow(x,[0.0,0.0],0.8);
JvDir_Handle=@(x,t)[zeros(size(x)) zeros(size(x))];

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
    h = CFL*mw(j);
    
    NewMesh = init_LEB(NewMesh);
    NewMesh=add_Edge2Elem(NewMesh);
    NewMesh=add_Patches(NewMesh);
    NewMesh = add_DGData(NewMesh);
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P3O3());  
    %M = assemMat_Mass1fD(NewMesh);
    
    % discrete pullback
    % Interpolation
    pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
    pbVp = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,h);
    
    % transported Mesh
    pMesh.Coordinates=pbVm(:,[1,2]);
    pMesh.Elements=Mesh.Elements;
    pMesh = add_Edges(pMesh);
    Loc = get_BdEdges(Mesh);
    pMesh.BdFlags = zeros(size(pMesh.Edges,1),1);
    pMesh.BdFlags(Loc) = -1;
    Mm = assemMat_W1F(pMesh,@MASS_W1F,MU_HANDLE, P3O3());

    P_i_dir = assemMat_SemiLag_W1F(NewMesh, pbVm);
    P_i_ad = assemMat_SemiLag_W1F(NewMesh, pbVp);
    
    AvJump=assemMat_Inn_LieW1F(NewMesh,@STIMA_Inn_AvgJumpW1F,Dir_Handle,gauleg(0,1,5));
    
    %[interpolation]
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
   
    %Dirichlet data
    [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),0);
    H_new=[H_Dir];
    %system matrix
    A = M;
    %righthand side
    Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),0);
    L=[Lf];
    L = h*M*L+M*P_i_dir* H_old-A*H_new;
    %solve system
    H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);

    %Dirichlet data
    [B_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),h);
    B_new=[B_Dir];
    %system matrix
    A = M;
    %righthand side
    Lg = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),0);
    L=Lg;
    %plot_Norm_W1F(Lg,NewMesh); colorbar;
    L = h*M*L+Mm*P_i_dir*B_old-A*B_new;
    %L = P_i_ad'*M* B_old-A*B_new;
    %solve system
    B_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
    %B_new = A\L;

   % plot_Norm_W1F(H_old,NewMesh);colorbar;
    plot_Norm_W1F((H_old-H_new)/h,NewMesh);colorbar;

    %plot_Norm_W1F(B_old,NewMesh);colorbar;
    plot_Norm_W1F((B_old-B_new)/h,NewMesh);colorbar;
   
    err(1,j)=L2Err_W1F_mod(NewMesh,(H_old-H_new)/h,P7O6(),@(x,varargin)zeros(size(x)),[0.0,0.0],0.5);
    err(2,j)=L2Err_W1F_mod(NewMesh,(B_old-B_new)/h,P7O6(),@(x,varargin)zeros(size(x)),[0.0,0.0],0.5)

    % patches 
    defMesh = NewMesh;
    defMesh.Coordinates = pbVp(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVp);
    P_p_dir = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    defMesh = NewMesh;
    defMesh.Coordinates = pbVm(:,[1 2]);
    intersec = aff_elems2(NewMesh, defMesh,pbVm);
    P_p_ad = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
    H_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    B_old = assemCochain_1f(NewMesh,H_Handle,gauleg(0,1,10),0);
    
    %Dirichlet data
    [H_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),0);
    H_new=[H_Dir];
    %system matrix
    A = M;
    %righthand side
    Lf = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),0);
    L=[Lf];
    L = h*M*L+P_p_dir* H_old-A*H_new;
    %solve system
    H_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);

    %Dirichlet data
    [B_Dir,FreeDofs] = assemDir_W1F(NewMesh,-1,H_Handle,gauleg(0,1,5),h);
    B_new=[B_Dir];
    %system matrix
    A = M;
    %righthand side
    Lg = assemCochain_1f(NewMesh, F_HANDLE,gauleg(0,1,10),0);
    L=Lg;
    %plot_Norm_W1F(Lg,NewMesh); colorbar;
    L = h*M*L+P_p_ad'* B_old-A*B_new;
    %L = P_i_ad'*M* B_old-A*B_new;
    %solve system
    B_new(FreeDofs,:) = A(FreeDofs,FreeDofs)\L(FreeDofs,:);
    %B_new = A\L;

   % plot_Norm_W1F(H_old,NewMesh);colorbar;
    plot_Norm_W1F((H_old-H_new)/h,NewMesh);colorbar;

    %plot_Norm_W1F(B_old,NewMesh);colorbar;
    plot_Norm_W1F((B_old-B_new)/h,NewMesh);colorbar;
    
    err(3,j)=L2Err_W1F_mod(NewMesh,(H_old-H_new)/h,P7O6(),@(x)zeros(size(x)),[0.0,0.0],0.5);
    err(4,j)=L2Err_W1F_mod(NewMesh,(B_old-B_new)/h,P7O6(),@(x)zeros(size(x)),[0.0,0.0],0.5)
end%  Mesh refinement

figure;
hold on;
plot(mw,err(1,:),'r-',mw,err(2,:),'g-',...
    mw,err(3,:),'r--',mw,err(4,:),'g--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf mesh width}','FontSize',12,'FontWeight','bold');
ylabel('{\bf L^2-error}','FontSize',12,'FontWeight','bold');
legend('direct Inter','adjoint Inter','direct Patch','adjoint Patch','Location','Northwest');