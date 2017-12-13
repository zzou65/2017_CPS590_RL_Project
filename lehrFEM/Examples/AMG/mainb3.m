function [] = mainb3()
% spectral radius of amg on locally refined meshes

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % initialize parameters
  
  theta0 = 0.5;
  itlvl0 = 1;
  cyc0 = 1;
  m0 = 1;
  smoother0 = @gs_smooth;
  tol0 = 0;
  maxit0 = 3;
  err0 = 'r';
  
  cyc1 = 1;
  m1 = 1;
  smoother1 = @gs_smooth;
  tol1 = 0;
  maxit1 = 1;
  
  numits = 30;
  
  total_iter = 9;
  initrefs = 12;
  
  % initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  [Ufull,Mesh] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,err0,theta0,itlvl0,cyc0,m0,smoother0,0,initrefs);
  
  % initialize data
  
  dofs = zeros(1,total_iter);
  rho_lmg = zeros(1,total_iter);
  rho_amg = zeros(1,total_iter);
  
  % refine mesh and calculate number of dofs and number of required
  % subspace corrections after each refinement
  
  for i=1:total_iter
    
    % refine mesh
    
    [Ufull,Mesh,flag,iter,ndof,err,timer,dofslvl,Afull,bfull,FDofs] = ...
      locmg_solve(Ufull,Mesh,@f_LShap,@g_D_LShap,err0,theta0,itlvl0,cyc0,m0,smoother0,tol0,maxit0);
    
    % determine number of degrees of freedom
    
    dofs(i) = ndof(end);
    
    % calculate convergence rate for multigrid using the power method;
    % the error propagation operator is E = I - BA
    
    Eerr = zeros(size(Ufull));
    Eerr(FDofs) = 1;
    u0 = zeros(size(Ufull));
    for j=1:numits
      err = Eerr;
      c = locmg_solve(u0,Mesh,Afull*err,Afull,u0,FDofs,[],cyc1,m1,smoother1,tol1,maxit1); % c = BAerr
      Eerr(FDofs) = err(FDofs) - c(FDofs);
    end
    rho_lmg(i) = (err(FDofs)'*Eerr(FDofs))/(err(FDofs)'*err(FDofs));
    
    % generate AMG data structure
    
    AMGOptions = AMGDefaultOptions;
    AMGOptions.mincoarse = 25;
    AMGOptions.pre.its = m1(1);
    AMGOptions.post.its = m1(end);
    
    A = Afull(FDofs,FDofs);

    AMGData = AMGSetup(A,AMGOptions);
    
    % calculate spectral radius of AMG with power method
    
    Eerr = ones(dofs(i),1);
    for j=1:numits
      err = Eerr;
      c = AMGVcycle(AMGData,A*err);
      Eerr = err - c;
    end
    rho_amg(i) = (err'*Eerr)/(err'*err);
    
  end

  % plot spectral radius vs number of degrees of freedom
  
  figure;
  plot(dofs,rho_lmg,'-^',dofs,rho_amg,'-+');
  hold on;
  set(gca,'XScale','log','YLim',[0 1]);
  grid('on');

  title('{\bf Dependence of Multigrid Convergence Rates on Refinement Level}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf spectral radius}')
  legend('Local Multigrid','Algebraic Multigrid','Location','NorthWest');

return
