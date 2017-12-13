% Run script for the maximal and minimal eigen values of stiffness and L2
% mass matrix on the quasi-uniform regularly refined mesh of unit square 
% domain using Lanzcos with Jacobian preconditioner.

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
    
    MAX_A = zeros(NREFS,1);
    MIN_A = zeros(NREFS,1);
    MAX_M = zeros(NREFS,1);
    MIN_M = zeros(NREFS,1);
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

    % Plot the solution
    
    plot(M_W,MIN_A,'-o',M_W,MAX_A,'-x')
    set(gca,'YScale','log','XScale','log')
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','NorthWest')
    title('\bf\fontsize{10}eigen values of Stiffness matrix on quasi-uniform refined mesh(using JAC PREC)')
    xlabel('\bf\fontsize{13}h_{min}')
    ylabel('\bf\fontsize{10}eigen value')
    grid on    
    print('-depsc','eigen_Qua_AJ.eps')
    
    figure
    plot(M_W,MIN_M,'-o',M_W,MAX_M,'-x')
    set(gca,'YScale','log','XScale','log')
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','NorthWest')
    title('\bf\fontsize{10}eigen values of MASS matrix on quasi-uniform refined mesh(using JAC PREC)')
    xlabel('\bf\fontsize{13}h_{min}')
    ylabel('\bf\fontsize{10}eigen value')
    set(gca,'YLim',[.1 5])
    grid on    
    print('-depsc','eigen_Qua_MJ.eps')
    
    % Clear memory
    
%     clear all
    
    
    