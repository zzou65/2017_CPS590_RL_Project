% Run script for the goal oriented error estimator.

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @(x,varargin)ones(size(x,1),1);   % Right hand side source term
  GD_HANDLE = @(x,varargin)zeros(size(x,1),1); % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap;                      % Exact siolution
  TOL = .2;                                   % Stopping criterion
  THETA = .3;                                   % Refinement strategy
  NITER = 5;                                    % Plotting levels
  PLOTTING = 0;                                 % plotting flag
  zero_fun = @(x,varargin) zeros(size(x,1),1);  % left hand side for dual problem

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
      
  iter = 0;
  Eta = TOL+1;
  while(Eta > TOL)
      
    iter = iter+1;
    nDofs(iter) = size(Mesh.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute the solution of the dual problem
    
    phi_interp = phi(Mesh.Coordinates);
    L_d = A'*phi_interp;
    [gF,FreeDofs]=assemDir_LFE(Mesh,-1,zero_fun);
    L_d = L_d-A*gF;
    A_T = A';
    gF(FreeDofs) = A_T(FreeDofs,FreeDofs)\L_d(FreeDofs);
    
    Mesh = add_Patches(Mesh);
    Eta_K = ErrEst_GOAL(gF,U,F_HANDLE,Mesh,P7O6(),QuadRule_1D);
    
    if PLOTTING
        plot_P0(Eta_K,Mesh);colorbar;
    end
    
    
    Eta(iter) = sqrt(sum(Eta_K.^2));
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    Mesh = refine_LEB(Mesh,MarkedElem);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta  :  %2.4e,  TOL  :  %2.4e\n',Eta(iter),TOL);
    
  end
  
  % Compute FE solution
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
  [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  L = L - A*U;
  if(~isempty(FreeDofs))
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
  end
    
  % Generate figures
  
  plot_LFE(U,Mesh);colorbar;
  plot_Mesh(Mesh,'as');
  
  fig = figure('Name','Goal oriented error estimator');
  plot(nDofs,Eta,'r-x');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Convergence rates for goal based error estimator}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  p = polyfit(log(nDofs),log(Eta),1);
  add_Slope(gca,'SouthWest',p(1));
    
  
  
  % Clear memory
  
  clear all;