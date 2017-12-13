% Run script for the residual based error estimator.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @f_LShap;          % Right hand side source term
  GD_HANDLE = @g_D_LShap;       % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap_H1S;  % Exact siolution
  TOL = .05;                    % Stopping criterion
  THETA = .5;                   % Refinement strategy
  NITER = 5;                    % Plotting levels 
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
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
    
    % Real and estimated errors
    
    Err(iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
    Mesh = add_Edge2Elem(Mesh);
    Eta_K = ErrEst_RES(U,Mesh,P7O6(),F_HANDLE);    
    Eta(iter) = sqrt(sum(Eta_K.^2));
       
    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
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
      
  % Generate figures
    
  fig = figure('Name','Residual based error estimator');
  plot(nDofs,Eta,'r-x', ...
       nDofs,Err,'b-x');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Convergence rates for residual based error estimator}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  legend('EST','Discr. err.','Location','NorthEast')
  p = polyfit(log(nDofs),log(Eta),1);
  add_Slope(gca,'SouthWest',p(1));
    
  fig = figure('Name','Residual based error estimator');
  plot(nDofs,Eta./Err,'r-x');
  grid('on');
  set(gca,'XScale','log');
  title('{\bf Ratio between estimated and discretization error}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Error ratio}');
    
  % Clear memory
  
  clear all;