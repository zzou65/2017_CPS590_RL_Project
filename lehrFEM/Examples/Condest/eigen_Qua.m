% Run script for the maximal and minimal eigen values of stiffness and L2
% mass matrix on the quasi-uniform regularly refined mesh of unit square domain

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NREFS = 6;
    Tol = 1e-3;
    Maxit_A = 10;
    Maxit_M = 10;
    
    % Preallocate memory
    
    MAX_A = zeros(NREFS,1);
    MIN_A = zeros(NREFS,1);
    MAX_M = zeros(NREFS,1);
    MIN_M = zeros(NREFS,1);
    M_W = zeros(NREFS,1);
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1; 1 -1; 1 1;-1 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh = refine_REG(Mesh);
  
    % Compute eigen values
    
    for iter = 1:NREFS
        
        Mesh = refine_REG(Mesh);
        Loc = get_BdEdges(Mesh);
        Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
        FixedPos = zeros(size(Mesh.Coordinates,1),1);
        FixedPos(Loc) = 1;
        JMesh = jiggle(Mesh,FixedPos);
        M_W(iter) = get_MeshMin(JMesh);
        A = assemMat_LFE(JMesh,@STIMA_Lapl_LFE);
        M = assemMat_LFE(JMesh,@MASS_LFE);
        nCoordinates = size(JMesh.Coordinates,1);
        Loc = get_BdEdges(JMesh);
        BDofs = unique(JMesh.Edges(Loc,:));
        FDofs = setdiff(1:nCoordinates,BDofs);
        A = A(FDofs,FDofs);     
                
        [dummy,MAX_A(iter)] = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'lm');
        [dummy,MIN_A(iter)] = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'sm');
        
        [dummy,MAX_M(iter)] = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'lm');
        [dummy,MIN_M(iter)] = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'sm');
                
    end

    % Plot the solution
    
    plot(M_W,MIN_A,'-o',M_W,MAX_A,'-x')
    set(gca,'YScale','log','XScale','log','XDir','reverse')
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','SouthWest')
    title('\bf\fontsize{10}eigen values of Stiffness matrix on quasi-uniform refined mesh')
    xlabel('\bf\fontsize{13}h_{min}')
    ylabel('\bf\fontsize{10}eigen value')
    grid on
    print('-depsc','eigen_Qua_A.eps')
    
    figure
    plot(M_W,MIN_M,'-o',M_W,MAX_M,'-x')
    set(gca,'YScale','log','XScale','log','XDir','reverse','YLim',[.1 10])
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','SouthWest')
    title('\bf\fontsize{10}eigen values of MASS matrix on quasi-uniform refined mesh')
    xlabel('\bf\fontsize{13}h_{min}')
    ylabel('\bf\fontsize{10}eigen value')
    grid on
    print('-depsc','eigen_Qua_M.eps')
    
    % Clear memory
    
%     clear all
    
    
    