% Run script for the condition number of stiffness and mass matrix with 
% bpx preconditioner on the unit square domain

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

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

    % Plot the solution
    
    plot(bpx_Cond_A,'-rx')
    hold on
    plot(bpx_Cond_M,'-b+')
    set(gca,'YScale','log')
    
   