% Run script for the maximal and minimal eigen values of stiffness and L2
% mass matrix on the distorted mesh of unit square domain using Lanzcos
% with Jacobian preconditioner.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NREFS = 3;
    NJIGGLE = 10;
    Tol = 1e-3;
    Maxit_A = 10;
    Maxit_M = 10;
    
    % Preallocate memory
    
    MAX_A = zeros(NREFS,1);
    MIN_A = zeros(NREFS,1);
    MAX_M = zeros(NREFS,1);
    MIN_M = zeros(NREFS,1);
    flag_A = ones(NREFS,1);
    flag_M = ones(NREFS,1);
    SR = zeros(NREFS,1);
    
    % Initialize mesh
    
    Mesh.Coordinates = [0 0; 1 0; 1 1;0 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;

    for i = 1:NREFS
        Mesh = refine_REG(Mesh);
    end


    Loc = get_BdEdges(Mesh);
    Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
    FixedPos = zeros(size(Mesh.Coordinates,1),1);
    FixedPos(Loc) = 1;
    nCoordinates = size(Mesh.Coordinates,1);
    BDofs = unique(Mesh.Edges(Loc,:));
    FDofs = setdiff(1:nCoordinates,BDofs);

    % Compute eigen values

    for iter = 1:NJIGGLE
        
        Mesh = jiggle(Mesh,FixedPos);
        SR(iter) = shape_reg(Mesh);
        A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
        M = assemMat_LFE(Mesh,@MASS_LFE);
        A = A(FDofs,FDofs);
        [MIN_A(iter) MAX_A(iter) flag_A(iter)] = eigen_Lanzcos(A,Tol,Maxit_A,@jac_prec,A);
        
        while( flag_A(iter) ~= 1 )
            Maxit_A = Maxit_A *2;
            [MIN_A(iter) MAX_A(iter) flag_A(iter)] = eigen_Lanzcos(A,Tol,Maxit_A,@jac_prec,A);
        end
        
        [MIN_M(iter) MAX_M(iter) flag_M(iter)] = eigen_Lanzcos(M,Tol,Maxit_M,@jac_prec,M);
        
        while( flag_M(iter) ~= 1 )
            Maxit_M = Maxit_M *2;
            [MIN_M(iter) MAX_M(iter) flag_M(iter)] = eigen_Lanzcos(M,Tol,Maxit_M,@jac_prec,M);
        end
        
        
    end

    % save the data table
    
    Table_A = [MAX_A MIN_A SR];    
    Table_M = [MAX_M MIN_M SR];
    save eigen_DisJ.mat Table_A Table_M
    
    % Clear memory
    
    clear all;
    