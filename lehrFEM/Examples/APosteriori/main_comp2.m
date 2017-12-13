% Run script for comparing the residual based error estimator, the goal
% oriented error estimator and uniform refinement

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  JIG=1;
  F_HANDLE = @f_LShap;                          % Right hand side source term
  GD_HANDLE = @g_D_LShap;                       % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap;                      % Exact solution
  UEX_Part = @(x,varargin)UEX_HANDLE(x).*indicator(x);% solution on part of the domain
  TOL = 10^-6;                                  % Stopping criterion
  MAX_DOFS = 10^5;                              % Stopping criterion
  THETA = .3;                                   % Refinement strategy
  NITER = 5;                                    % Plotting levels
  PLOTTING = 0;                                 % plotting flag
  NREFS = 4;                                    % number of red refinements
  zero_fun = @(x,varargin) zeros(size(x,1),1);  % left hand side for dual problem

  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = refine_REG(Mesh);
  Mesh = init_LEB(Mesh);
  
  Mesh_REG = Mesh;
  Mesh_GOAL = Mesh;
  Mesh_RES = Mesh;
  
  
  fprintf('\nStandard Mesh refinement \n\n');
  
  iter = 0;
  nDofs_REG = 0;
  while nDofs_REG < MAX_DOFS
    iter = iter+1;
    
    nDofs_REG(iter) = size(Mesh_REG.Coordinates,1);  
    
    % Mesh preprocessing
    
    switch(JIG)
    case 1
      NewMesh_REG = Mesh_REG;      
    case 2
      Loc = get_BdEdges(Mesh_REG);
      Loc = unique([Mesh_REG.Edges(Loc,1); Mesh_REG.Edges(Loc,2)]);
      FixedPos = zeros(size(Mesh_REG.Coordinates,1),1);
      FixedPos(Loc) = 1;
      NewMesh_REG = jiggle(Mesh_REG,FixedPos);   
    case 3
      Loc = get_BdEdges(Mesh_REG);
      Loc = unique([Mesh_REG.Edges(Loc,1); Mesh_REG.Edges(Loc,2)]);
      FixedPos = zeros(size(Mesh_REG.Coordinates,1),1);
      FixedPos(Loc) = 1;
      NewMesh_REG = smooth(Mesh_REG,FixedPos);
    end
    Mesh_REG = NewMesh_REG;
      
    % solve boundary value problem
    A = assemMat_LFE(Mesh_REG,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh_REG,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh_REG,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute nodal interpolator of the fem solution on [0.5 1]^2
    
    U_part = U.*indicator(Mesh_REG.Coordinates);
    
    % compute L1 error on the square of interest
    
    err_REG(iter) = comp2_err(Mesh_REG,U_part,P7O6(),UEX_Part);
    err_REG_L1(iter) = L1Err_LFE(Mesh_REG,U_part,P7O6(),UEX_Part);
    
    % perform uniform mesh refinement
    
    tmp_Mesh = Mesh_REG;
    Mesh_REG = refine_REG(Mesh_REG);
    
  end
  
  Mesh_REG = tmp_Mesh;
  
  
  
  fprintf('\nResidual based adaptive mesh refinement \n\n')
  
  iter = 0;
  nDofs_RES = 0;
  while(nDofs_RES < MAX_DOFS)
  
    iter = iter+1;
    nDofs_RES(iter) = size(Mesh_RES.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh_RES,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh_RES,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh_RES,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute nodal interpolator of the fem solution on [0.5 1]^2
    
    U_part = U.*indicator(Mesh_RES.Coordinates);
    
    % compute L1 error on the square of interest
    
    err_RES(iter) = comp2_err(Mesh_RES,U_part,P7O6(),UEX_Part);
    err_RES_L1(iter) = L1Err_LFE(Mesh_RES,U_part,P7O6(),UEX_Part);
    
    % Real and estimated errors
    
    Mesh_RES = add_Edge2Elem(Mesh_RES);
    Eta_res_K = ErrEst_RES(U,Mesh_RES,P7O6(),F_HANDLE);    
    Eta_res(iter) = sqrt(sum(Eta_res_K.^2));
       
    % Mark elements for refinement
    
    Eta_max = max(Eta_res_K);
    MarkedElem = find(Eta_res_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh_RES = Mesh_RES;
    Mesh_RES = refine_LEB(Mesh_RES,MarkedElem);
    
    % Update mesh data structure
    
    Mesh_RES = add_Edges(Mesh_RES);
    Mesh_RES = rmfield(Mesh_RES,'BdFlags');
    Loc = get_BdEdges(Mesh_RES);
    Mesh_RES.BdFlags = zeros(size(Mesh_RES.Edges,1),1); 
    Mesh_RES.BdFlags(Loc) = -1; 
        
    % Print out information
    
    fprintf('Eta  :  %2.4e,  Dofs  :  %2.4e\n',Eta_res(iter),nDofs_RES(iter));
  
  end
  
    
  
  
  
  fprintf('\nGoal oriented adaptive mesh refinement \n\n');
  
  QuadRule_1D = gauleg(0,1,5);
      
  iter = 0;
  nDofs_GOAL = 0;
  while(nDofs_GOAL < MAX_DOFS)
      
    iter = iter+1;
    nDofs_GOAL(iter) = size(Mesh_GOAL.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh_GOAL,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh_GOAL,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh_GOAL,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute nodal interpolator of the fem solution on [0.5 1]^2
    
    U_part = U.*indicator(Mesh_GOAL.Coordinates);
    
    % compute L1 error on the square of interest
    
    err_GOAL(iter) = comp2_err(Mesh_GOAL,U_part,P7O6(),UEX_Part);
    err_GOAL_L1(iter) = L1Err_LFE(Mesh_GOAL,U_part,P7O6(),UEX_Part);
    
    % compute the solution of the dual problem
    
    L_d = assemLoad_LFE(Mesh_GOAL,P7O6(),@indicator);
    [gF,FreeDofs]=assemDir_LFE(Mesh_GOAL,-1,zero_fun);
    L_d = L_d-A*gF;
    A_T = A';
    gF(FreeDofs) = A_T(FreeDofs,FreeDofs)\L_d(FreeDofs);
    
    Mesh_GOAL = add_Patches(Mesh_GOAL);
    Eta_goal_K = ErrEst_GOAL(gF,U,F_HANDLE,Mesh_GOAL,P7O6(),QuadRule_1D);
    
    
    Eta_goal(iter) = sqrt(sum(Eta_goal_K.^2));
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_goal_K);
    MarkedElem = find(Eta_goal_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh_goal = Mesh_GOAL;
    Mesh_GOAL = refine_LEB(Mesh_GOAL,MarkedElem);
    
    % Update mesh data structure
    
    Mesh_GOAL = add_Edges(Mesh_GOAL);
    Mesh_GOAL = rmfield(Mesh_GOAL,'BdFlags');
    Loc = get_BdEdges(Mesh_GOAL);
    Mesh_GOAL.BdFlags = zeros(size(Mesh_GOAL.Edges,1),1); 
    Mesh_GOAL.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta  :  %2.4e,  Dofs  :  %2.4e\n',Eta_goal(iter),nDofs_GOAL(iter));
    
  end
  
  
    
  % Generate figures
  
  plot_Mesh(Mesh_REG,'as');
  FileName = 'Uniform.eps';
  print('-depsc2', FileName); 
  
  plot_Mesh(Mesh_RES,'as');
  FileName = 'Residual.eps';
  print('-depsc2', FileName); 
  
  plot_Mesh(Mesh_GOAL,'as');
  FileName = 'Goal.eps';
  print('-depsc2', FileName); 
  
  
  fig = figure('Name','Discretization errors');
  plot(nDofs_REG,abs(err_REG),'r-x',...
      nDofs_RES,abs(err_RES),'-bx',...
      nDofs_GOAL,abs(err_GOAL),'-gx');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Convergence rates for monitored functional}');
  legend('Uniform ref.','Residual ref.','Goal ref.','Location','NorthEast')
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  p = polyfit(log(nDofs_REG),log(abs(err_REG)),1);
  add_Slope(gca,'SouthWest',p(1),'r-');
  p = polyfit(log(nDofs_RES),log(abs(err_RES)),1);
  add_Slope(gca,'East',p(1),'b-');
  p = polyfit(log(nDofs_GOAL),log(abs(err_GOAL)),1);
  add_Slope(gca,'North',p(1),'g-');
  
  FileName = 'Functional.eps';
  print('-depsc2', FileName); 
  
    
  fig = figure('Name','L^1 Discretization errors');
  plot(nDofs_REG,abs(err_REG_L1),'r-x',...
      nDofs_RES,abs(err_RES_L1),'-bx',...
      nDofs_GOAL,abs(err_GOAL_L1),'-gx');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Convergence rates for L^1 norm on square of interest}');
  legend('Uniform ref.','Residual ref.','Goal ref.','Location','NorthEast')
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  p = polyfit(log(nDofs_REG),log(abs(err_REG_L1)),1);
  add_Slope(gca,'SouthWest',p(1),'r-');
  p = polyfit(log(nDofs_RES),log(abs(err_RES_L1)),1);
  add_Slope(gca,'East',p(1),'b-');
  p = polyfit(log(nDofs_GOAL),log(abs(err_GOAL_L1)),1);
  add_Slope(gca,'North',p(1),'g-');
  
  FileName = 'L^1_norm.eps';
  print('-depsc2', FileName); 
  
  % Clear memory
  
  clear all;