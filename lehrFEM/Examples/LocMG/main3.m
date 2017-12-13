function [] = main3()
% distribution of degrees of freedom for refinement towards a circle

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % choose error estimators
 
  reconstruct = false;
  point = true;

  % initialize parameters
  
  theta = 0.5;
  itlvl = 1;
  cyc = 1;
  m = 1;
  smoother = @gs_smooth;
  tol = 0;
  maxit = 10;
  
  its = [1,4,7,10];
  
  % initialize mesh

  Mesh = load_Mesh('Coord_Sqr_Big.dat','Elem_Sqr_Big.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:2
    Mesh = refine_REG(Mesh);
  end
  
  % initialize rhs
  
  assemLoad = @(mesh) assemLoad_LFE_S1(mesh,gauleg(0,1,2));
  g_D = @(x,varargin) zeros(size(x,1),1);
  
  % do reconstruction-based error estimator stuff
  
  if(reconstruct)
    
    % run multigrid with reconstruction-based error estimator

    [U_c,Mesh_c,flag_c,iter_c,ndof_c,err_c,timer_c,dofslvl_c] = locmg_solve([],Mesh,{assemLoad},g_D,'c',theta,itlvl,cyc,m,smoother,tol,maxit);

    LVL_c = max(Mesh_c.VertLevel);

    % plot distribution of degrees of freedom for reconstruction-based error
    % estimator

    figure;
    plot((0:LVL_c)',dofslvl_c(:,its));
    xlabel('\bf level');
    ylabel('\bf number of degrees of freedom');
    title('\bf Distribution of Degrees of Freedom for Reconstruction-Based Error Estimator');
    set(gca,'XLim',[0,LVL_c],'XTick',(0:LVL_c));
    grid on;
    legend('iteration 1','iteration 4','iteration 7','iteration 10','Location','NorthWest');
    
  end % reconstruction-based error estimator
  
  
  % do point-refinement
  
  if(point)
    
    % run multigrid with hierarchical error estimator

    [U_p,Mesh_p,flag_p,iter_p,ndof_p,err_p,timer_p,dofslvl_p] = locmg_solve([],Mesh,{assemLoad},g_D,@circref,theta,itlvl,cyc,m,smoother,tol,maxit);

    LVL_p = max(Mesh_p.VertLevel);

    % plot distribution of degrees of freedom for residual-based error
    % estimator

    figure;
    plot((0:LVL_p)',dofslvl_p(:,its));
    xlabel('\bf level');
    ylabel('\bf number of degrees of freedom');
    title('\bf Distribution of Degrees of Freedom for Refinement Towards a Circle');
    set(gca,'XLim',[0,LVL_p],'XTick',(0:LVL_p));
    grid on;
    legend('iteration 1','iteration 4','iteration 7','iteration 10','Location','NorthWest');
    
    plot_Mesh(Mesh_p,'as');
    
  end % reconstruction-based error estimator
  
  