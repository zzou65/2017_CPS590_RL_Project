function [] = mainb()
% convergence of amg on locally refined meshes

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize parameters
  
  refs = 0:25:100;          % number of local refinements
  drefs = diff([0,refs]);
  num = length(refs);
  
  tol = 1e-6;
  maxit = 50;

  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:5
    Mesh = refine_REG(Mesh);
  end
  Mesh = init_LEB(Mesh);
  
  % Initialize solution
  
  iter = nan(1,num);
  resvec = cell(1,num);
  
  for j=1:num
  
    % Do local refinements towards 0

    for i=1:drefs(j)
      Mesh = add_Patches(Mesh);
      eta = pointref([],Mesh);
      marked = find(eta);
      Mesh = refine_LEB(Mesh,marked);
      Mesh = add_Edges(Mesh);
    end
    
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;

    % Assemble stiffness matrix

    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    b = assemLoad_LFE(Mesh,P7O6(),@f_LShap);
    [U,FDofs] = assemDir_LFE(Mesh,-1,@g_D_LShap);
    b = b - A*U;
    
    b = b(FDofs);
    A = A(FDofs,FDofs);
    
    % Generate AMG data structure

    AMGOptions = AMGDefaultOptions;
    AMGOptions.mincoarse = 10;
    
    AMGData = AMGSetup(A,AMGOptions);
    
    % Run AMG solver
    
    [x,flag,relres,iter(j),resvec{j}] = amg_solve(zeros(size(b)),b,[],AMGData,tol,maxit);
    
  end % for loop with index j
  
  % Plot residuals against iteration
  
  figure;
  lgd = cell(1,num);
  for j=1:num
    semilogy(1:iter(j),resvec{j});
    hold all;
    lgd{j} = sprintf('%.0f refinements',refs(j));
  end
  legend(lgd{:},'Location','NorthEast');
  title('\bf Convergence of AMG on Locally Refined Meshes');
  xlabel('\bf iteration');
  ylabel('\bf residual');
  grid on;
  
return