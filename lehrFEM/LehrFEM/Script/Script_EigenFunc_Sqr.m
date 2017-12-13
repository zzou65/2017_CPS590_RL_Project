% Run script for the Error of eigenfunction in the unit square domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    GUEX_HANDLE_1 = @(x,varargin)[-pi/2*sin(pi/2*x(:,1)).*cos(pi/2*x(:,2)),...
                                  -pi/2*cos(pi/2*x(:,1)).*sin(pi/2*x(:,2))];
    GUEX_HANDLE_2 = @(x,varargin)[pi*cos(pi*x(:,1)).*sin(pi*x(:,2))...
                                  pi*sin(pi*x(:,1)).*cos(pi*x(:,2))];
    NEigen = 4;                                     % Number of the eigenvalue
    NREFS = 7;                                      % Number of red refinement steps
    Select = [1 4];

    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_Func = zeros(NREFS,2);
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh);
    
    for iter = 1:NREFS
        
        % Do REG refine
        
        Mesh = refine_REG(Mesh);
        
        % Get mesh width
        
        M_W(iter) = get_MeshWidth(Mesh);
        
        % Assemble stiffness matrix and mass matrix

        A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
        M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());


        % Incorporate Dirichlet boundary data

        [U,FreeNodes] = assemDir_LFE(Mesh,-1,GD_HANDLE);
        A = A(FreeNodes,FreeNodes);
        M = M(FreeNodes,FreeNodes);

        % Solve eigenvalue problem

        [ V d ] = eigs(A,M,NEigen,'sm');
        
        for i = 1:NEigen

            U(FreeNodes,i) = V(:,i);
          
        end

        % Scaling the solution to fit the exact one

        U(:,NEigen) = abs(U(:,NEigen));
        [dummy pos] = max(U(:,NEigen-3));
        vert = Mesh.Coordinates(pos,:);
        if (vert(1)*vert(2)<0)
            U(:,NEigen-3) = -U(:,NEigen-3);
        end

        % Compute the error

        Err_Func(iter,1) = H1SErr_LFE(Mesh,U(:,NEigen),P7O6(),GUEX_HANDLE_1);
        Err_Func(iter,2) = H1SErr_LFE(Mesh,U(:,NEigen-3),P7O6(),GUEX_HANDLE_2);
        
    end
    
    figure;
    plot(M_W,abs(Err_Func(:,1)),'-x',...
         M_W,abs(Err_Func(:,2)),'-^'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for eigenfunction errors(Square) LFE');
    xlabel('\bfh');
    ylabel('\bfenergy error(H1S)');
    legend('EigenFunc.1','EigenFunc.4',...
            'Location','NorthEast')
    p = polyfit(log(M_W),log(abs(Err_Func(:,1))),1);
    add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_EigFunc_sqr.eps');
    !gv rate_EigFunc_sqr.eps &
    
    % Clear memory
    
    clear all;
 