% Run script for Crouzeix-Raviart finite element solver.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland   

    % Initialize constant
  % Clear memory
    
    clear all;
    
    NREFS = 4;
    F_Handle = @(x,varargin)2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
    GD_Handle = @(x,varargin)cos(pi*x(:,1)).*cos(pi*x(:,2));
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1; 1 1;-1 1];
    Mesh.Elements =[1 2 4;2 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;
    
    for i=1:NREFS
       
        Mesh = refine_REG(Mesh);
        
    end
    
    % Assemble stiffness matrix and load vector
    
    A = assemMat_CR(Mesh,@STIMA_Lapl_CR);
    L = assemLoad_CR(Mesh,P7O6(),F_Handle);
    
    % Incorporate Dirichlet boundary data

    [U,FreeDofs] = assemDir_CR(Mesh,-1,GD_Handle);
    L = L - A*U;

    % Solve the linear system
 
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % Plot the solution
    
    plot_CR(U,Mesh);
    colorbar;
    
  
    
    
    