% Run script for the eigenvalue problem in the unit disk domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    NEigen = 10;                    % Number of the eigenvalue
    Select = [1 2 4 6 7 9 11 13 15 16 18 20 22 24];
    ZBessel = [2.4048 3.8317 5.1356 6.3802 7.5883 8.7715 5.5201...
              7.0156 8.4172 9.7610 11.0647 8.6537 10.1735 9.93611];
    ZBessel = sort(ZBessel)';
    L = 1:14;
    Lambda = ZBessel.^2;
    
    % Initialize mesh
    
    NREFS = 3;             % Number of inital mesh refinements
    DHANDLE = @dist_circ;  % Signed distance function
    C = [0 0];             % Center of the circle
    R = 1;                 % Radius of the circle

    % Initialize mesh

    Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    for i = 1:NREFS
        Mesh = refine_REG(Mesh,DHANDLE,C,R);
    end
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_ParBd(Mesh,DHANDLE,C,R);

    % Assemble mass matrix and load vector

    A = assemMat_PBD(Mesh,@STIMA_Lapl_PBD,P7O6());
    M = assemMat_PBD(Mesh,@MASS_PBD,P7O6());

    % Incorporate Dirichlet boundary conditions

    [U,FreeNodes] = assemDir_PBD(Mesh,-1,GD_HANDLE);
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
        
        if (NEigen+1-i > 0)
            
        plot_PBD(U(:,i),Mesh);
        colorbar;
        title(['{\bf Solution:' int2str(NEigen+1-i) '  eigenvalue = ' num2str(d(i,i)^.5) '}']);
        end
    end
    
%     % Output the .eps file
%     
%     d = diag(d);
%     d = d(NEigen:-1:1);
%     LambdaH = d(Select(1:6));
%     plot(L',log(Lambda)-5.8,'-',...
%          L',log((LambdaH-Lambda)./Lambda),'^')
%     grid on
%     title('{\bf Errors in the higher eigenvalues of a disk}')
%     xlabel('{\bf Number of eigenvalue}')
%     ylabel('{\bf rel. Error [log]}')
%     legend(['\bflog\lambda_{l}+c'],['\bflog((\lambda^{h}_{l}-\lambda_{l})/\lambda_{l})'],...
%             'Location','East')
    
    % Clear memory
    
%     clear all;