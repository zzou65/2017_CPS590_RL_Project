% Run script for the Error of eigenvalue problem in the L-shaped domain

%   Copyright 2008-2008 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    NEigen = 6;                                     % Number of the eigenvalue
    NP = 12;                                        % Number of spectral refinement steps
    P0 = 2;
    Select = [2 3 4 6];
    Lambda = [1.47562182408 3.53403136678 9.86960440109 11.3894793979];
    
    % Preallocate memory
    P = 2+(1:NP);
    D_LambdaH = zeros(NP,4);
    
    % Initialize mesh
    Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    Mesh = add_Edge2Elem(Mesh);
    
    % Loop over polynomial degrees
    for j=1:length(P)
      p = P(j);
      
      % Determine dof data
      EDofs = (p-1)*ones(size(Mesh.Edges,1),1);
      CDofs = (p-1)*(p-2)/2*ones(size(Mesh.Elements,1),1);
      Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);
      QuadRule = Duffy(TProd(gauleg(0,1,2*p)));
      Shap = shap_hp(QuadRule.x,p);
      
      % Assemble mass and stiffness matrices
      A = assemMat_hp(Mesh,Elem2Dof,@STIMA_Lapl_hp,QuadRule,Shap);
      M = assemMat_hp(Mesh,Elem2Dof,@MASS_hp,QuadRule,Shap);
      
      % Compute eigenvalues
%       [V,d] = eigs(A,M,NEigen,'sm');
      [V,d] = eig(full(A),full(M));
      d = sort(diag(d));
      
      % Determine errors in eigenvalues
      D_LambdaH(j,:) = d(Select).' - Lambda;
      
    end
    
    figure;
    plot(P,abs(D_LambdaH(:,1)),'-x',...
         P,abs(D_LambdaH(:,2)),'-^',...
         P,abs(D_LambdaH(:,3)),'-h',...
         P,abs(D_LambdaH(:,4)),'-o'); 
    grid on
    set(gca,'YScale','log','XTick',P,'XLim',[min(P),max(P)],'YLim',[1e-12,1e0]);
    title('\bfConvergence rate for eigenvalue errors (L-shape) p-FEM');
    xlabel('\bf p');
    ylabel('\bf|\lambda^{p}_{l}-\lambda_{l}|');
    legend('EigenVal.1','EigenVal.2','EigenVal.3','EigenVal.4',...
            'Location','SouthWest')
    print('-depsc', 'Eig_LShap_pFEM.eps');
    !gv Eig_LShap_pFEM.eps & 
    
    figure;
    plot(P,abs(D_LambdaH(:,1)),'-x',...
         P,abs(D_LambdaH(:,2)),'-^',...
         P,abs(D_LambdaH(:,3)),'-h',...
         P,abs(D_LambdaH(:,4)),'-o'); 
    grid on
    set(gca,'XScale','log','YScale','log','XTick',P,'XLim',[min(P),max(P)],'YLim',[1e-6,1e0]);
    title('\bfConvergence rate for eigenvalue errors (L-shape) p-FEM');
    xlabel('\bf p');
    ylabel('\bf|\lambda^{p}_{l}-\lambda_{l}|');
    legend('EigenVal.1','EigenVal.2','EigenVal.3','EigenVal.4',...
            'Location','NorthEast')
    p = polyfit(log(P'),log(abs(D_LambdaH(:,2))),1);
    add_Slope(gca,'SouthWest',p(1));
    p1 = polyfit(log(P'),log(abs(D_LambdaH(:,1))),1);
    add_Slope(gca,'North',p1(1));
    print('-depsc', 'Eig_LShap_pFEM_log.eps');
    !gv Eig_LShap_pFEM_log.eps & 
    
    % Clear memory
    clear all;   