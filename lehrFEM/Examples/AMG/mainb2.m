function [] = mainb2()
% convergence of amg on locally refined meshes

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Initialize parameters
  
  m = [1,1];
  maxit = 20;
  tol = 1e-6;

  % Initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);

  % Generate mesh (with incorporated mulrigrid structure) with local multigrid
  
  [Ufull,Mesh,flag,iter,ndof,err,timer,dofslvl,Afull,bfull,FDofs] = ...
    locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'r',0.5,1,1,m,@gs_smooth,0,21);
  
  % Incorporate Dirichlet boundary conditions and calculate
  
  U = Ufull;
  U(FDofs) = 0;
  b = bfull - Afull*U;
  b = b(FDofs);
  A = Afull(FDofs,FDofs);
  U_ex = A\b;
  
  err0 = sqrt(U_ex'*A*U_ex);
  
  %% Run (unnested) AMG on locally refined mesh
  
  % Generate AMG data structure

  AMGOptions = AMGDefaultOptions;
  AMGOptions.mincoarse = 25;
  AMGOptions.pre.its = m(1);
  AMGOptions.post.its = m(end);

  AMGData = AMGSetup(A,AMGOptions);
    
  % Initialize AMG
    
  errors_amg = nan(1,maxit);
  U_amg = zeros(size(b));
  
  % Run AMG solver
  
  for i=1:maxit
    U_amg = amg_solve(U_amg,b,[],AMGData,0,1);
    errors_amg(i) = sqrt((U_amg-U_ex)'*A*(U_amg-U_ex));
    if(errors_amg(i) < tol)
      break;
    end
  end
  iter_amg = i;
  errors_amg = [err0,errors_amg(1:iter_amg)];
  
  %% Run (unnested) multigrid solver on locally refined mesh
  
  % Initialize

  Ufull_ex = U;
  Ufull_ex(FDofs) = U_ex;
  
  % Run multigrid solver (locmg_solve code)
  
  [U_mg,Mesh,flag,iter_mg,ndof,errors_mg] = ...
    locmg_solve([],Mesh,bfull,Afull,U,FDofs,Ufull_ex,1,m,@gs_smooth,tol,maxit);

  errors_mg = [err0,errors_mg];
  
  %% Count 1-d subspace corrections in V-cycles
  
  % for locmg_solve
  
  dofs_mg = dofslvl(:,end);
  numssc_mg = (m(1)+m(end))*sum(dofs_mg);
  
  % for amg
  
  L = length(AMGData);
  dofs_amg = zeros(1,L);
  for i=1:L
    dofs_amg(i) = size(AMGData{L-i+1}.A,1);
  end
  numssc_amg = (m(1)+m(end))*sum(dofs_amg);
  
  %% Plot convergence
  
  figure;
  semilogy(numssc_mg*(0:iter_mg),errors_mg,'-^',...
    numssc_amg*(0:iter_amg),errors_amg,'-+');
  legend('Local Multigrid','Algebraic Multigrid','Location','NorthEast');
  title('\bf Convergence of Multigrid Methods on Locally Refined Meshes');
  xlabel('\bf number of subspace corrections');
  ylabel('\bf error (energy norm)');
  grid on;
  
  fprintf('dofs : %d\n',ndof(end));
  fprintf('LVL  : %d\n',max(Mesh.VertLevel));

return