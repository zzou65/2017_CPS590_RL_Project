function [] = main2()
% distribution of degrees of freedom for refinement towards a point

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize parameters
  
  theta = 0.5;
  itlvl = 1;
  cyc = 1;
  m = 1;
  smoother = @gs_smooth;
  tol = 0;
  maxit = 20;
  
  its = [1,6,11,15,18,20];

  % initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:2
    Mesh = refine_REG(Mesh);
  end
    
  % run multigrid with residual-based error estimator

  [U_r,Mesh_r,flag_r,iter_r,ndof_r,err_r,timer_r,dofslvl_r] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'r',theta,itlvl,cyc,m,smoother,tol,maxit);

  LVL_r = max(Mesh_r.VertLevel);

  % plot distribution of degrees of freedom for residual-based error
  % estimator

  figure;
  plot((0:LVL_r)',dofslvl_r(:,its));
  xlabel('\bf level');
  ylabel('\bf number of degrees of freedom');
  title('\bf Distribution of Degrees of Freedom for Residual-Based Error Estimator');
  set(gca,'XLim',[0,LVL_r]);
  grid on;
  legend('iteration 1','iteration 6','iteration 11','iteration 15','iteration 18','iteration 20','Location','NorthEast');

  % plot initial mesh and refines mesh
  
  plot_Mesh(Mesh,'as');
  
  plot_Mesh(Mesh_r,'as');
  
  
 

  
return
