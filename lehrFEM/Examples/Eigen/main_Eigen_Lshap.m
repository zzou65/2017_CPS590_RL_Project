% Run script for the eigenvalue problem in the L-Shaped domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant

    NEigen = 3;                          % Number of the eigenvalue
    NREFS = 5;                            % Number of red refinement steps

    % Initialize mesh
    
    Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    for i = 1:NREFS
        Mesh = refine_REG(Mesh);
    end   
    plot_Mesh(Mesh,'sa')
    
    % Assemble stiffness matrix and mass matrix
    
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
    M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    
    % Solve eigenvalue problem
    
    [ V d ] = eigs(A,M,NEigen,'sm');
    
    for i = 1:NEigen
        
        U(:,i) = V(:,i);
        norm(A*V(:,i)-d(i,i)*M*V(:,i))
      
    end
      
    % Plot the solution
    
    for i = NEigen:-1:1
        
        plot_LFE(U(:,i),Mesh);
        colorbar;
        title(['{\bf Solution:' int2str(NEigen+1-i) '  eigenvalue = ' num2str(d(i,i)) '}']);
        
    end
    
    % Clear memory
    
    clear all;
 
    
    