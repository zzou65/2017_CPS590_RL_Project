% Run script for the maximal and minimal eigen values of stiffness matrix
% preconditioned by the Mass matrix on the distorted mesh of unit square 
% domain.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NREFS = 4;
    NJIGGLE = 10;
    Tol = 1e-3;
    Maxit_A = 10;
    Maxit_M = 10;
    
    % Preallocate memory
    
    MAX_A = zeros(NREFS,1);
    MIN_A = zeros(NREFS,1);
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
        M = M(FDofs,FDofs);
        
        [dummy,MIN_A(iter)] = eigs(A,M,1,'sm');
        [dummy,MAX_A(iter)] = eigs(A,M,1,'lm');        
        
    end

    % save the data table
    
    Table_A = [MAX_A MIN_A SR];
    save stima_eig_Dis.mat Table_A
    
    % Clear memory
    
    clear all;
    