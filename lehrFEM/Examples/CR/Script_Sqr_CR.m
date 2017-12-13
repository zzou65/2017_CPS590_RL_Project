% Run script for energy error convergence rate for CR finite element on the
% unit square domain.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland   

    % Initialize constant
    
    NREFS = 7;
    F_Handle = @(x,varargin)2*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2));
    GD_Handle = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));
    UEX_Handle = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));
    UEX_GD_Handle = @(x,varargin)pi*[cos(pi*x(:,1)).*sin(pi*x(:,2)) ...
                                     sin(pi*x(:,1)).*cos(pi*x(:,2))];
    
    % Preallocate memory
    
    Err_L2_CR = zeros(NREFS,1);
    Err_H1S_CR = zeros(NREFS,1);
    Err_L2_LFE = zeros(NREFS,1);
    Err_H1S_LFE = zeros(NREFS,1);
    M_W = zeros(NREFS,1);
    
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1; 1 1;-1 1];
    Mesh.Elements =[1 2 4;2 3 4];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;
    Mesh = refine_REG(Mesh);
    
    for i = 1:NREFS
       
        Mesh = refine_REG(Mesh);
        M_W(i) = get_MeshWidth(Mesh);
        
        % Compute CR solution
        
        A_CR = assemMat_CR(Mesh,@STIMA_Lapl_CR);
        L_CR = assemLoad_CR(Mesh,P7O6(),F_Handle);
        [U_CR,FreeDofs] = assemDir_CR(Mesh,-1,GD_Handle);
        L_CR = L_CR - A_CR*U_CR;
        U_CR(FreeDofs) = A_CR(FreeDofs,FreeDofs)\L_CR(FreeDofs);
        
        % Compute error
        
        Err_L2_CR(i) = L2Err_CR(Mesh,U_CR,P7O6(),UEX_Handle);
        Err_H1S_CR(i) = H1SErr_CR(Mesh,U_CR,P7O6(),UEX_GD_Handle);
        
        % Compute LFE solution
        
        A_LFE = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
        L_LFE = assemLoad_LFE(Mesh,P7O6(),F_Handle);
        [U_LFE,FreeDofs] = assemDir_LFE(Mesh,-1,GD_Handle);
        L_LFE = L_LFE - A_LFE*U_LFE;
        U_LFE(FreeDofs) = A_LFE(FreeDofs,FreeDofs)\L_LFE(FreeDofs);
        
        % Compute error
        
        Err_L2_LFE(i) = L2Err_LFE(Mesh,U_LFE,P7O6(),UEX_Handle);
        Err_H1S_LFE(i) = H1SErr_LFE(Mesh,U_LFE,P7O6(),UEX_GD_Handle);
        
    end
    
   
    figure;
    plot(M_W,Err_L2_CR,'-^r',M_W,Err_H1S_CR,'-xb',...
         M_W,Err_L2_LFE,'-s',M_W,Err_H1S_LFE,'-o'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for L2& Broken H1S error on unit square');
    xlabel('\bf\fontsize{14}h');
    ylabel('\bf\fontsize{14}Energy error');
    legend('CR: L^{2}-error','CR: H^{1}-error',...
           'PL: L^{2}-error','PL: H^{1}-error'); 
    p = polyfit(log(M_W),log(Err_H1S_CR),1);
    add_Slope(gca,'West',p(1));
    p = polyfit(log(M_W),log(Err_L2_CR),1);
    add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_Sqr_CR.eps');
    !gv rate_Sqr_CR.eps &
    
    % Clear memory
    
    clear all
    
    