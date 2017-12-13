function [] = main4()
% convergence of local multigrid on fixed refined meshes

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
  maxit0 = 2;
  err0 = 'r';
  
  cyc1 = 1;
  m1 = 1;
  smoother1 = @gs_smooth;
  tol1_c = 0.1;
  maxit1 = 100;
  
  total_iter = 23;
  initrefs = 2;
  
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
  numssc = zeros(1,total_iter);
  
  % refine mesh and calculate number of dofs and number of required
  % subspace corrections after each refinement
  
  for i=1:total_iter
    
    % refine mesh
    
    [Ufull,Mesh,flag,iter,ndof,err,timer,dofslvl,Afull,bfull,FDofs] = ...
      locmg_solve(Ufull,Mesh,@f_LShap,@g_D_LShap,err0,theta0,itlvl0,cyc0,m0,smoother0,tol0,maxit0);
    
    % determine number of degrees of freedom and cost of V-cycle
    
    dofs(i) = ndof(end);
    numssc(i) = (m1(1)+m1(end))*sum(dofslvl(:,end));
    
    % incorporate Dirichlet boundary conditions
    
    U_bd = Ufull;
    U_bd(FDofs) = 0;
    b = bfull - Afull*U_bd;
    b = b(FDofs);
    A = Afull(FDofs,FDofs);
    U_ex = A\b;
    
    U_ex_full = U_bd;
    U_ex_full(FDofs) = U_ex;
    
    % calculate discretization error
    
    disc_err = H1SErr_LFE(Mesh,U_ex_full,P7O6(),@grad_uex_LShap);
    tol1 = tol1_c*disc_err;
    
    % run multigrid solver
    
    [U_mg,Mesh,flag,iter] = locmg_solve([],Mesh,bfull,Afull,U_bd,FDofs,U_ex_full,cyc1,m1,smoother1,tol1,maxit1);
    
    numssc(i) = numssc(i)*flag*iter;
    
  end

  % plot number of subspace corrections vs number of degrees of freedom
  
  figure;
  plot(dofs,numssc,'o');
  hold on;
  set(gca,'XScale','log','YScale','log');
  grid('on');

  title('{\bf Complexity of Local Mutligrid}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf one-dimensonal subspace corrections}')

  p = polyfit(log(dofs),log(numssc),1);
  x = get(gca,'XLim');
  y = get(gca,'YLim');
  plot(x,exp(polyval(p,log(x))));
  set(gca,'YLim',y);
  add_Slope(gca,'NorthWest',p(1));


return