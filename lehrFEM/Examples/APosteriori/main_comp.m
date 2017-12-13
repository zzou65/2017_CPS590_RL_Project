% Compare residual based and goal oriented error estimators for obtaining
% the outer boundary flux

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @(x,varargin)ones(size(x,1),1);   % Right hand side source term
  GD_HANDLE = @(x,varargin) zeros(size(x,1),1); % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap;                      % Exact siolution
  TOL = .09;                                    % Stopping criterion
  THETA_res = .5; THETA_goal = .1;             % Refinement strategy
  NITER = 5;                                    % Plotting levels
  PLOTTING = 0;                                 % plotting flag
  zero_fun = @(x,varargin) zeros(size(x,1),1);  % left hand side for dual problem
  limit = load('limit.mat');                    % load approximate limit of output functional
  limit = limit.limit;
  MAX_NDOFS = 10^5;
  coeff = 1;

  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr2.dat','Elem_Sqr2.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = refine_REG(Mesh);
  Mesh = init_LEB(Mesh);
  
  QuadRule_1D = gauleg(0,1,5);
  
  Mesh_goal = Mesh;
  Mesh_res = Mesh;
      
  
  fprintf('\nResidual based adaptive mesh refinement \n\n')
  
  iter = 0;
  nDofs_res = 0;
  while(nDofs_res < MAX_NDOFS)
      
    iter = iter+1;
    
    nDofs_res(iter) = size(Mesh_res.Coordinates,1);  
    
    % compute finite element solution for residual based estimator
    
    A = assemMat_LFE(Mesh_res,@STIMA_Lapl_LFE);
    B = assemMat_LFE(Mesh_res,@MASS_LFE);    
    L = assemLoad_LFE(Mesh_res,P7O6(),F_HANDLE);
    [U_res,FreeDofs] = assemDir_LFE(Mesh_res,-1,GD_HANDLE);
    L = L - A*U_res;
    if(~isempty(FreeDofs))
      U_res(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute error estimator
    
    Mesh_res = add_Edge2Elem(Mesh_res);
    Eta_K_res = ErrEst_RES(U_res,Mesh_res,P7O6(),F_HANDLE);    
    Eta_res(iter) = sqrt(sum(Eta_K_res.^2));
    
    % compute finite element approximation of cutoff function phi
    
    [phi_interp,FreeDofs] = assemDir_LFE(Mesh_res,-1,@gd_phi);
    M = A+coeff*B;
    L = -M*phi_interp;
    phi_interp(FreeDofs) = M(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % compute error output functional
    
    F_res(iter) = U_res'*A*phi_interp - integr(Mesh_res,phi_interp,P7O6);
    Err_res(iter) = (F_res(iter));
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_K_res);
    MarkedElem_res = find(Eta_K_res > Eta_max*THETA_res);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh_res = Mesh_res;
    Mesh_res = refine_LEB(Mesh_res,MarkedElem_res);
    
    % Update mesh data structure
    
    Mesh_res = add_Edges(Mesh_res);
    Mesh_res = rmfield(Mesh_res,'BdFlags');
    Loc = get_BdEdges(Mesh_res);
    Mesh_res.BdFlags = zeros(size(Mesh_res.Edges,1),1); 
    Mesh_res.BdFlags(Loc) = -1;  
    
    % Print out information
    
    fprintf('Dofs  :  %2.4e,  F(u-u_n)  :  %2.4e\n',nDofs_res(iter),abs(Err_res(iter)-limit));
  
  end
  
  
  fprintf('\nGoal oriented adaptive mesh refinement \n\n');  
  
  iter = 0;
  nDofs_goal = 0;
  while(nDofs_goal < MAX_NDOFS)
    iter = iter+1;
    nDofs_goal(iter) = size(Mesh_goal.Coordinates,1);  
    
     % Compute FE solution for the goal oriented mesh
  
    A = assemMat_LFE(Mesh_goal,@STIMA_Lapl_LFE);
    B = assemMat_LFE(Mesh_goal,@MASS_LFE);    
    L = assemLoad_LFE(Mesh_goal,P7O6(),F_HANDLE);
    [U_goal,FreeDofs] = assemDir_LFE(Mesh_goal,-1,GD_HANDLE);
    L = L - A*U_goal;
    if(~isempty(FreeDofs))
      U_goal(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute finite element approximation of cutoff function phi
    
    [phi_interp,FreeDofs] = assemDir_LFE(Mesh_goal,-1,@gd_phi);
    M = A+coeff*B;
    L = -M*phi_interp;
    phi_interp(FreeDofs) = M(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % compute error output functional
    
    F_goal(iter) = U_goal'*A*phi_interp - integr(Mesh_goal,phi_interp,P7O6);
    Err_goal(iter) = (F_goal(iter));
    
    % compute the solution of the dual problem
    
    L_d = A'*phi_interp;
    [gF,FreeDofs]=assemDir_LFE(Mesh_goal,-1,zero_fun);
    L_d = L_d-A*gF;
    A_T = A';
    gF(FreeDofs) = A_T(FreeDofs,FreeDofs)\L_d(FreeDofs);
    
    % cpmpute error estimator
    
    Mesh_goal = add_Patches(Mesh_goal);
    Eta_K_goal = ErrEst_GOAL(gF,U_goal,F_HANDLE,Mesh_goal,P7O6(),QuadRule_1D);
    Eta_goal(iter) = sqrt(sum(Eta_K_goal.^2));
    
        
    % Mark elements for refinement
    
    Eta_max = max(Eta_K_goal);
    MarkedElem_goal = find(Eta_K_goal > Eta_max*THETA_goal);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh_goal = Mesh_goal;
    Mesh_goal = refine_LEB(Mesh_goal,MarkedElem_goal);
    
    % Update mesh data structure
    
    Mesh_goal = add_Edges(Mesh_goal);
    Mesh_goal = rmfield(Mesh_goal,'BdFlags');
    Loc = get_BdEdges(Mesh_goal);
    Mesh_goal.BdFlags = zeros(size(Mesh_goal.Edges,1),1); 
    Mesh_goal.BdFlags(Loc) = -1; 
    
    fprintf('Dofs  :  %2.4e,  F(u-u_n)  :  %2.4e\n',nDofs_goal(iter),abs(Err_goal(iter)-limit));
  end
  
  Err_goal = abs(Err_goal - limit);
  Err_res = abs(Err_res - limit);
  
  
  % Generate figures
  
  plot_LFE(U_goal,tmp_Mesh_goal);colorbar;
  plot_Mesh(tmp_Mesh_goal);
  FileName = 'Mesh_goal.eps';
  print('-depsc2', FileName); 
  
  plot_LFE(U_res,tmp_Mesh_res);colorbar;
  plot_Mesh(tmp_Mesh_res);
  FileName = 'Mesh_res.eps';
  print('-depsc2', FileName); 
  
  fig = figure('Name','Comparison of error estimators');
  plot(nDofs_goal,Err_goal,'r-x',...
      nDofs_res,Err_res,'-bx');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Comparison of error estimators}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  legend('Goal oriented','Residual based','Location','NorthEast')
  p = polyfit(log(nDofs_goal),log(Err_goal),1);
  add_Slope(gca,'SouthWest',p(1),'-r');
  p = polyfit(log(nDofs_res),log(Err_res),1);
  add_Slope(gca,'West',p(1),'-b');
  
  FileName = 'compare.eps';
  print('-depsc2', FileName); 
    
  
  
  % Clear memory
  
  clear;