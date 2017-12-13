% Run script for the condition number of stiffness and mass matrix on the
% unit square domain

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NREFS = 7;
    Tol = 1e-3;
    Maxit_A = 10;
    Maxit_M = 10;
    
    % Preallocate memory
    
    Cond_A = zeros(NREFS,1);
    Cond_M = zeros(NREFS,1);
    flag_A = ones(NREFS,1);
    flag_M = ones(NREFS,1);
    M_W = zeros(NREFS,1);
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1; 1 -1; 1 1;-1 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
  
    % Compute condition number
    
    for iter = 1:NREFS
        
        Mesh = refine_REG(Mesh);
        M_W(iter) = get_MeshWidth(Mesh);
        A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
        M = assemMat_LFE(Mesh,@MASS_LFE);
        nCoordinates = size(Mesh.Coordinates,1);
        Loc = get_BdEdges(Mesh);
        BDofs = unique(Mesh.Edges(Loc,:));
        FDofs = setdiff(1:nCoordinates,BDofs);
        A = A(FDofs,FDofs);
        [Cond_A(iter) flag_A(iter)] = condest_Lanzcos(A,Tol,Maxit_A);
        
        while( flag_A(iter) ~= 1 )
            Maxit_A = Maxit_A *2;
            [Cond_A(iter) flag_A(iter)] = condest_Lanzcos(A,Tol,Maxit_A);
        end
        
        [Cond_M(iter) flag_M(iter)] = condest_Lanzcos(M,Tol,Maxit_M);
        
        while( flag_M(iter) ~= 1 )
            Maxit_M = Maxit_M *2;
            [Cond_M(iter) flag_M(iter)] = condest_Lanzcos(M,Tol,Maxit_M);
        end
        
        
    end
    
    save Cond.mat M_W Cond_A Cond_M
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % With preconditioner 
    
    % Initialize constants
    
    clear
    IREFS = 2;
    NREFS = 6;
    Tol = 1e-3;
    Maxit = 10;
    
    % Preallocate memory
    
    bpx_Cond_A = zeros(NREFS+1,1);
    bpx_Cond_M = zeros(NREFS+1,1);
    
    % Initialize mesh
    
    CMesh.Coordinates = [-1 -1; 1 -1; 1 1;-1 1];
    CMesh.Elements = [1 2 3;1 3 4];
    CMesh.ElemFlag = ones(size(CMesh.Elements,1),1);
    CMesh = add_Edges(CMesh);
    Loc = get_BdEdges(CMesh);
    CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
    CMesh.BdFlags(Loc) = -1;

    for i = 1:IREFS
        CMesh = refine_REG(CMesh);
    end
    
    % Compute non-Dirichlet vertices

    Loc = get_BdEdges(CMesh);
    DDofs = unique([CMesh.Edges(Loc,1); CMesh.Edges(Loc,2)]);
    CFDofs = setdiff(1:size(CMesh.Coordinates,1),DDofs);

    % Generate multigrid data structure

    for i = 1:NREFS

        % Compute stiffness matrix and load vector

        A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
        A = A(CFDofs,CFDofs);
        MG_Data_S.D{i} = diag(A);
        M = assemMat_LFE(CMesh,@MASS_LFE);
        MG_Data_M.D{i} = diag(M);
        
        if i>1
            flag = 1;

            [bpx_Cond_A(i) flag] = condest_Lanzcos(A,Tol,Maxit,@bpx_prec,MG_Data_S);

            while( flag ~= 1 )
                Maxit = Maxit *2
                [bpx_Cond_A(i) flag] = condest_Lanzcos(A,Tol,Maxit,@bpx_prec,MG_Data_S);
            end
        end
        
        if i>1
            flag = 1;

            [bpx_Cond_M(i) flag] = condest_Lanzcos(M,Tol,Maxit,@bpx_prec,MG_Data_M);

            while( flag ~= 1 )
                Maxit = Maxit *2
                [bpx_Cond_M(i) flag] = condest_Lanzcos(M,Tol,Maxit,@bpx_prec,MG_Data_M);
            end
        end
        
        % Refine the mesh and compute prolongation matrix

        FMesh = refine_REG(CMesh);
        P = get_PMat_LFE(CMesh,FMesh);
        MG_Data_M.P{i} = P;
        Loc = get_BdEdges(FMesh);
        DDofs = unique(FMesh.Edges(Loc,:));
        FFDofs = setdiff(1:size(FMesh.Coordinates,1),DDofs);
        MG_Data_S.P{i} = P(FFDofs,CFDofs);

        % Update coarse mesh

        CMesh = FMesh;
        CFDofs = FFDofs;

    end

    % on finest mesh
    
    A = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
    Loc = get_BdEdges(CMesh);
    DDofs = unique([CMesh.Edges(Loc,1); CMesh.Edges(Loc,2)]);
    CFDofs = setdiff(1:size(CMesh.Coordinates,1),DDofs);
    A = A(CFDofs,CFDofs);
    MG_Data_S.D{NREFS+1} = diag(A);
    
    M = assemMat_LFE(CMesh,@MASS_LFE);
    MG_Data_M.D{NREFS+1} = diag(M);
    
    % Compute condition number

    flag = 1;

    [bpx_Cond_A(NREFS+1) flag] = condest_Lanzcos(A,Tol,Maxit,@bpx_prec,MG_Data_S);

    while( flag ~= 1 )
        Maxit = Maxit *2
        [bpx_Cond_A(NREFS+1) flag] = condest_Lanzcos(A,Tol,Maxit,@bpx_prec,MG_Data_S);
    end

    flag = 1;

    [bpx_Cond_M(NREFS+1) flag] = condest_Lanzcos(M,Tol,Maxit,@bpx_prec,MG_Data_M);

    while( flag ~= 1 )
        Maxit = Maxit *2
        [bpx_Cond_M(NREFS+1) flag] = condest_Lanzcos(M,Tol,Maxit,@bpx_prec,MG_Data_M);
    end
    
    save bpx_Cond.mat bpx_Cond_A bpx_Cond_M
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Generate figures
    
    % Load data
    
    load Cond.mat
    load bpx_Cond.mat
    
    % Generate figure
    
    plot(M_W,Cond_A,'-*',M_W,Cond_M,'-x',M_W,bpx_Cond_A,'-o',M_W,bpx_Cond_M,'->')
    
    set(gca,'XScale','log','YScale','log','XDir','reverse')
    
    legend('\bf\fontsize{10}Cond S','\bf\fontsize{10}Cond M',...
           '\bf\fontsize{10}BPX S','\bf\fontsize{10}BPX M','location','Northwest');
    xlabel('\bf\fontsize{12}Mesh width')
    ylabel('\bf\fontsize{12}condition number')
    title('\bf\fontsize{10}Condition number of Stiffness and Mass matrix')
    
    % Output .eps file
    
    print('-depsc','rate_Cond.eps')
    
    % Clear memory
    
    delete *.mat
    clear all