% Run script for weak regularization W1F finite element solver.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 2;
    U_Handle = @(x,varargin)ones(size(x,1),1);
%     F_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))];
%     GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    F_Handle = @(x,varargin)[zeros(size(x,1),1) zeros(size(x,1),1)];
    GD_Handle = @(x,varargin)[-x(:,2) x(:,1)];
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=-1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
    
        Mesh = refine_REG(Mesh);
        
    end
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    
    tmp = [];
    nCoordinates = size(Mesh.Coordinates,1);
    [IC,JC,C] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
    B = assemMat_WRegW1F(Mesh,@STIMA_WReg_W1F);
    D = assemMat_LFE(Mesh,@MASS_Lump_LFE);
    
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == -1);
    DNodes = unique([Mesh.Edges(DEdges,1); Mesh.Edges(DEdges,2)]);
    tmp = [tmp; DNodes];
    FreeDofs = setdiff(1:nCoordinates,tmp);
    
    B = B(:,FreeDofs);
    D = D(FreeDofs,FreeDofs);
    
    T = B*inv(D)*transpose(B);
    [IT,JT,T] = find(T);
    A = sparse([IC;IT],[JC;JT],[C;T]);
    L = assemLoad_W1F(Mesh,P7O6(),F_Handle);

    % Incorporate Dirichlet boundary data

    [U,FreeDofs] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
    L = L - A*U;
    
    % Solve the system
    
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % Plot the solution
    
    plot_W1F(U,Mesh);
    
    % Clear memory
    
%     clear all
    