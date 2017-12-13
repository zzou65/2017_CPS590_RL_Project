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

    % Plot the solution
    
    plot(1:NREFS,Cond_A,'-o',1:NREFS,Cond_M,'-x')
    set(gca,'YScale','log')
    legend('\bf\fontsize{10}cond_A','\bf\fontsize{10}cond_M')
    title('\bf\fontsize{10}precondition number')
    xlabel('\bf\fontsize{10}Step of refinement')
    ylabel('\bf\fontsize{10}Precondition number')
    