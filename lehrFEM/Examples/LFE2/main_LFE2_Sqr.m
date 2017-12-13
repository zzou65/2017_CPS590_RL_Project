% Run script for LFE2 finite element solver.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 5;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))]+[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=-1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
    
        Mesh = refine_REG(Mesh);
        
    end
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    
    [IC,JC,C] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
    [IM,JM,M] = assemMat_LFE2(Mesh,@MASS_LFE2);
    A = sparse([IC;IM],[JC;JM],[C;M]);
    L = assemLoad_LFE2(Mesh,P7O6(),F_Handle);
    
    % Incorporate Dirichlet boundary data
    
    [U,FreeDofs] = assemDir_LFE2(Mesh,-1,GD_Handle);
    L = L - A*U;
    
    % Solve the system
    
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % Plot the solution
    
    plot_DomBd(Mesh,'as')
    hold on
    plot_LFE2(U,Mesh);
    
    % Clear memory
    
%     clear all
    
    