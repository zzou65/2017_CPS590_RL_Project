% Run script for L2 error convergence rate for CR finite element solver.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland   

    % Initialize constant
    
    NREFS = 7;
    F_Handle = @(x,varargin)2*pi^2*cos(pi*x(:,1)).*cos(pi*x(:,2));
    GD_Handle = @(x,varargin)cos(pi*x(:,1)).*cos(pi*x(:,2));
    UEX_Handle = @(x,varargin)cos(pi*x(:,1)).*cos(pi*x(:,2));
    % Preallocate memory
    
    Err_L2 = zeros(NREFS,1);
    M_W = zeros(NREFS,1);
    
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
        M_W(i) = get_MeshWidth(Mesh);
        
        % Compute numerical solution
        
        A = assemMat_CR(Mesh,@STIMA_Lapl_CR);
        L = assemLoad_CR(Mesh,P7O6(),F_Handle);
        [U,FreeDofs] = assemDir_CR(Mesh,-1,GD_Handle);
        L = L - A*U;
        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
        % Compute error
        
        Err_L2(i) = L2Err_CR(Mesh,U,P7O6(),UEX_Handle);
        
    end
    
    % Plot the convergence rate
    
    figure;
    plot(M_W,Err_L2,'-xb'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for L2 error on unit square(CR)');
    xlabel('\bfh');
    ylabel('\bfL2 error');
    p = polyfit(log(M_W),log(Err_L2),1);
    add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_Sqr_L2_CR.eps');
    !gv rate_Sqr_L2_CR.eps &
    
    % Clear memory
    
    clear all;
    
    