% Run script for the maximal and minimal eigen values of stiffness and L2
% mass matrix on the locally refined mesh of unit square domain

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NREFS = 5;
    Tol = 1e-3;
    Maxit_A = 10;
    Maxit_M = 10;
    
    % Preallocate memory
    
    MAX_A = zeros(NREFS,1);
    MIN_A = zeros(NREFS,1);
    MAX_M = zeros(NREFS,1);
    MIN_M = zeros(NREFS,1);
    M_W = zeros(NREFS,1);
    
    % Initialize mesh
   
    Mesh.Coordinates = [0 0; 1 0; 1 1;0 1];
    Mesh.Elements = [1 2 3;1 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh = refine_REG(Mesh);
    Mesh = init_LEB(Mesh);
    
    % Compute eigen values
    
    for iter = 1:NREFS
        
        [Index,dummy1,dummy2] = find(Mesh.Elements == 1);
        Index = unique(Index);
        Mesh = refine_LEB(Mesh,Index);
	
	% added by R.H.: plotting of mesh
	figure('name',sprintf('mesh on level %d',iter));
	plot_Mesh(Mesh,'as');
	print('-depsc',sprintf('eigen_Loc_mesh%d',iter));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
        M_W(iter) = get_MeshMin(Mesh);
        A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
        M = assemMat_LFE(Mesh,@MASS_LFE);
        nCoordinates = size(Mesh.Coordinates,1);
        Loc = get_BdEdges(Mesh);
        BDofs = unique(Mesh.Edges(Loc,:));
        FDofs = setdiff(1:nCoordinates,BDofs);
        A = A(FDofs,FDofs);
                
%        [dummy,MAX_A(iter)] = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'lm');
%        [dummy,MIN_A(iter)] = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'sm');
        
        MIN_A(iter) = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'sm');
        MAX_A(iter) = eigs(A,spdiags(diag(A),0,size(A,1),size(A,2)),1,'lm');

%        [dummy,MAX_M(iter)] = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'lm');
%        [dummy,MIN_M(iter)] = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'sm');

        MIN_M(iter) = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'sm');
        MAX_M(iter) = eigs(M,spdiags(diag(M),0,size(M,1),size(M,2)),1,'lm');
    end

    % Plot the solution
    
    plot(M_W,MIN_A,'-o',M_W,MAX_A,'-x')
    set(gca,'YScale','log','XScale','log','XDir','reverse')
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','SouthWest')
    title('\bf\fontsize{10}eigen values of scaled stiffness matrix (locally refined mesh)')
    xlabel('\bf\fontsize{13}h_{min}')
    ylabel('\bf\fontsize{10}eigen value')
    grid on
    print('-depsc','eigen_Loc_A.eps')
        
    figure
    plot(M_W,MIN_M,'-o',M_W,MAX_M,'-x')
    set(gca,'YScale','log','XScale','log','XDir','reverse');
    legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','SouthWest');
    title('\bf\fontsize{10}eigen values of scaled mas matrix (locally refined mesh)');
    xlabel('\bf\fontsize{13}h_{min}');
    ylabel('\bf\fontsize{10}eigen value');
    axis([0.9*min(M_W) 1.1*max(M_W) 0.1 4]);
% =======
%     set(gca,'YScale','log','XScale','log','XDir','reverse','YLim',[.1 10])
%     legend('\bf\fontsize{10}min','\bf\fontsize{10}max','location','SouthWest')
%     title('\bf\fontsize{10}eigen values of MASS matrix on locally refined mesh')
%     xlabel('\bf\fontsize{13}h_{min}')
%     ylabel('\bf\fontsize{10}eigen value')
% >>>>>>> .r717
    grid on
    print('-depsc','eigen_Loc_M.eps')
    
    % Clear memory
    
    % clear all
    
    
    