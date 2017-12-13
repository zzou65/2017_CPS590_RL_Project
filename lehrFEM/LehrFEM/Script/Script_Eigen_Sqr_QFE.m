% Run script for the Error of eigenvalue problem in the unit square domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    NEigen = 6;                                     % Number of the eigenvalue
    NREFS = 6;                                      % Number of red refinement steps
    Select = [1 2 4 5];
    Lambda = [2 5 8 10]*pi^2/4;
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    D_LambdaH = zeros(NREFS,4);
    
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

        A = assemMat_QFE(Mesh,@STIMA_Lapl_QFE,P7O6());
        M = assemMat_QFE(Mesh,@MASS_QFE,P7O6());


        % Incorporate Dirichlet boundary data

        [U,FreeNodes] = assemDir_QFE(Mesh,-1,GD_HANDLE);
        A = A(FreeNodes,FreeNodes);
        M = M(FreeNodes,FreeNodes);

        % Solve eigenvalue problem

        [ V d ] = eigs(A,M,NEigen,'sm');
        
        for i = 1:NEigen

            U(FreeNodes,i) = V(:,i);
            norm(A*V(:,i)-d(i,i)*M*V(:,i))

        end

        % Plot one of the solution

        if( iter == NREFS )
            
            plot_QFE(U(:,NEigen),Mesh);
            colorbar;
            title(['{\bf Solution: 1  '                          ...  
                   'eigenvalue =' num2str(d(NEigen,NEigen),4)    ... 
                   ', meshwidth = ' num2str(M_W(iter),3) '}' ]);
            print('-depsc', 'Eig_func1_sqr_QFE.eps');
            !gv Eig_func1_sqr_QFE.eps &            
            
            plot_QFE(U(:,NEigen-1),Mesh);
            colorbar;
            title(['{\bf Solution: 2  '                           ...  
                   'eigenvalue =' num2str(d(NEigen-1,NEigen-1),4) ... 
                   ', meshwidth = ' num2str(M_W(iter),3) '}' ]);
            print('-depsc', 'Eig_func2_sqr_QFE.eps');
            !gv Eig_func2_sqr_QFE.eps &   
            
        end

        % Compute the error

        d = diag(d);
        d = d(6:-1:1);
        D_LambdaH(iter,:) = d(Select)'-Lambda;
    
    end
    
    figure;
    plot(M_W,abs(D_LambdaH(:,1)),'-x',...
         M_W,abs(D_LambdaH(:,2)),'-^',...
         M_W,abs(D_LambdaH(:,3)),'-h',...
         M_W,abs(D_LambdaH(:,4)),'-o'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for eigenvalue errors(Square) QFE');
    xlabel('\bfh');
    ylabel('\bf|\lambda^{h}_{l}-\lambda_{l}|');
    legend('EigenVal.1','EigenVal.2','EigenVal.3','EigenVal.4',...
            'Location','NorthEast')
    p = polyfit(log(M_W),log(abs(D_LambdaH(:,4))),1);
    add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_Eig_sqr_QFE.eps');
    !gv rate_Eig_sqr_QFE.eps &
    
    % Clear memory
    
    clear all;
 