% Run script for the eigenvalue problem in the unit square domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    NEigen = 2;                                  % Number of the eigenvalue
    NREFS = 6;                                    % Number of red refinement steps

    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    for i = 1:NREFS
        Mesh = refine_REG(Mesh);
    end   
    
    % Assemble stiffness matrix and mass matrix
    
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
    M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    
    
    % Incorporate Neumann boundary data
    
    [U,FreeNodes] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    A = A(FreeNodes,FreeNodes);
    M = M(FreeNodes,FreeNodes);

    % Solve eigenvalue problem
    
    [ V d ] = eigs(A,M,NEigen,'sm');
    
    for i = 1:NEigen
        
        U(FreeNodes,i) = V(:,i);
        norm(A*V(:,i)-d(i,i)*M*V(:,i))
      
    end
      
    % Plot the solution
    
    for i = NEigen:-1:1
        
        plot_LFE(U(:,i),Mesh);
        colorbar;
        title(['{\bf Solution:' int2str(NEigen+1-i) '  eigenvalue = ' num2str(d(i,i)) '}']);
        
    end
    
    % Clear memory
    
%     clear all;
 
    
    