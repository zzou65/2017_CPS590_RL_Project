% test difference of mass matrix and P matrix for small time steps and
% different computational methods for P

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


clear Mesh;
%clear all;
NREFS_init = 0;     % Number of uniform red refinements
NREFS =2;
JIG =1;

MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin) 1;

CFL =2.^-[3:7]; %

v1=0.1;
v2=1;
DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
DIR2=@(x)v2*sin(pi.*x(:,2)).*sin(pi.*x(:,1));
D1DIR2=@(x)v2*pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
D2DIR2=@(x)v2*pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));
DIR2=@(x)v2*(1-x(:,2).^2).*(1-x(:,1).^2);
D1DIR2=@(x)v2*(1-x(:,2).^2).*(-2*x(:,1));
D2DIR2=@(x)v2*(-2*x(:,2)).*(1-x(:,1).^2);
% DIR2=@(x)v2*(1-x(:,1));
% D1DIR2=@(x)-v2*ones(size(x,1),1);
% D2DIR2=@(x)0*(x(:,2));

Dir_Handle=@(x,t)[DIR1(x), DIR2(x)];
JvDir_Handle=@(x,t)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
JvDirT_Handle=@(x,t)[D1DIR1(x) D2DIR1(x) D1DIR2(x) D2DIR2(x)];

% Load mesh from file
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');% initial mesh

% Add edge data structure
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

%Mesh=refine_REG(Mesh);
% refine
for i=1:NREFS
    
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

    %ID1= assemMat_W1F(NewMesh,@STIMA_ContrRot,Dir_Handle, P7O6());
    %DI1 = assemMat_W1F(NewMesh,@STIMA_GradContr,Dir_Handle, gauleg(0,1,3));
    LieVol = assemMat_W1F(NewMesh,@STIMA_Lie,Dir_Handle,JvDir_Handle);
    LieInn = assemMat_Lie_Inn_W1F(NewMesh,Dir_Handle);
%    LieInn2 = assemMat_Lie_Inn_W1F2(NewMesh,Dir_Handle);

     SP = SparsePattern(NewMesh);
     SP=1-(SP>0);
    
    for j=1:length(CFL)

        mw = get_MeshWidth(NewMesh);
        h = CFL(j)*mw/norm([v1,v2]);
        directions = Dir_Handle(NewMesh.Coordinates,h);
        
        % Quadrature
        %pbB = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,-h);
        %P_qbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);
        %pbB = trace_bcenters(NewMesh, Dir_Handle,JvDir_Handle,h);
        %P_q_adbary = assemMat_SemiLagQuad_W1F_bary(NewMesh, pbB);

        %pbVm = trace_vertices_W(NewMesh,Dir_Handle,JvDir_Handle,-h);
        
        %plot_transportMesh(NewMesh,pbVm,'east');
        %plot_MeshPatch(Mesh,pbVm,'');
        
        %P_i = assemMat_SemiLag_W1F(NewMesh, pbVm);

        M = assemMat_W1F(NewMesh,@MASS_W1F,MU_HANDLE, P1O2());

        pbV = trace_vertices(NewMesh,h*directions);
        NewMesh = init_LEB(NewMesh);
        NewMesh = add_Patches(NewMesh);
        NewMesh = add_Edge2Elem(NewMesh);
        defMesh = NewMesh;
        defMesh.Coordinates = pbV(:,[1 2]);
        intersec = aff_elems2(NewMesh, defMesh,pbV);

        %test2(NewMesh,defMesh,intersec);
        
        P_q = assemMat_SemiLagQuad_W1F_patches(NewMesh, defMesh,intersec);
        P_qVol = assemMat_SemiLagQuad_W1F_patchesVol(NewMesh, defMesh,intersec);
        P_qInn = assemMat_SemiLagQuad_W1F_patchesInn(NewMesh, defMesh,intersec);
        M_q = assemMat_MASSQuad_W1F_patches(NewMesh, defMesh,intersec);

        %pbV = trace_vertices(NewMesh,-h*directions);
        %defMesh.Coordinates = pbV(:,[1 2]);
        %intersec = aff_elems2(NewMesh, defMesh,pbV);
        %P_q_ad = assemMat_SemiLagQuad_W1F_patches(NewMesh,
        %defMesh,intersec)
        
        %C=SP.*(M_q-P_q)/h;
        %SLieInn=SP.*LieInn;
        
        set(gcf, 'renderer', 'ZBuffer');
        % err(i,j,:)=[abs(C(5,4)-SLieInn(5,4))];
        
        err(i,j,:)=[norm(full(M_q-P_q)/h) norm(full(M_q-P_qVol)/h) norm(full(-P_qInn)/h)];
        err_q(i,j,:)=[norm(full((M_q-P_qVol)/h-LieVol)) norm(full(-P_qInn/h-LieInn))];
        %err_qbary(i,j,:)=[norm(full(M-P_qbary)) norm(full(M-P_qbary)/h) norm(full(M-P_qbary)/h^2)];
        Clim=[min([min(min((M_q-P_q)/h)),min(min(LieInn+LieVol))]) max([max(max((M_q-P_q)/h)),max(max(LieInn+LieVol))])];
        figure;
        subplot(2,2,1); imagesc(((M_q-P_qVol)/h-LieVol));set(gca,'CLIM',Clim);colorbar;view([0,90]); title('(M_q-P_q)/h-LieVol');       
        subplot(2,2,2); imagesc((LieInn));set(gca,'CLIM',Clim);colorbar;view([0,90]); title('Lie_{Inn}');
        Clim=[min([min(min((-P_qInn)/h)),min(min(LieInn))]) max([max(max((-P_qInn)/h)),max(max(LieInn))])];
        subplot(2,2,3); imagesc(-P_qInn/h);set(gca,'CLIM',Clim);colorbar;view([0,90]);title('(P_qInn)');
        subplot(2,2,4); imagesc(LieInn);set(gca,'CLIM',Clim);colorbar;view([0,90]); title('Lie_{Inn}');
        
        %figure;
        %subplot(1,2,1); imagesc((M_q-P_qVol-P_qInn/h));colorbar;view([0,90]);title('(M_q-P_q)');
        %subplot(1,2,2); imagesc(LieInn+LieVol);colorbar;view([0,90]); title('Lie_{Inn}');
        
    end
    Mesh=refine_REG(Mesh);
end

err_q