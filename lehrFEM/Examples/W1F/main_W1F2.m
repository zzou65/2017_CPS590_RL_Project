% run script for W1F finite element solver including a convective term

%   Copyright 2009-2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 4;
    MU_HANDLE=@(x,varargin)1;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)(pi^2+1)*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    W_Handle = @(x,varargin)[zeros(size(x,1),1) 1*ones(size(x,1),1)];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
%     F_Handle = @(x,varargin)[2-x(:,2) x(:,1)];
%     GD_Handle = @(x,varargin)[-x(:,2) x(:,1)];

    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
    
        Mesh = refine_REG(Mesh);
        
    end
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    t = cputime;
    [IC,JC,C] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
    [IM,JM,M] = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P3O3());
    [ICo,JCo,Co] = assemMat_W1F(Mesh,@CONV_Curl_W1F,W_Handle, P3O3());
    A = sparse([IC;IM;ICo],[JC;JM;JCo],[C;M;Co]);
    L = assemLoad_W1F(Mesh,P7O6(),F_Handle);
    
    % Incorporate Dirichlet boundary data
    
    [U,FreeDofs] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
    L = L - A*U;
    
    % Solve the system
    
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    fprintf('Runtime of direct solver [s]  :  %f\n',cputime-t);
    % Plot the solution
    
%     plot_Mesh(Mesh,'as')
%     hold on
    plot_W1F(U,Mesh);
    %plot_Norm_W1F(U,Mesh);

    
    % Clear memory
    
     clear all