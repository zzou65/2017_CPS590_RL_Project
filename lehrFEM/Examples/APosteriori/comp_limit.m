  % script for computing the limit of the monitored functional F as
  % accurate as possible

  
  %   Copyright 2009 Christoph Wiesmeyr
  %   SAM - Seminar for Applied Mathematics
  %   ETH-Zentrum
  %   CH-8092 Zurich, Switzerland

  F_HANDLE = @(x,varargin)ones(size(x,1),1);    % Right hand side source term
  GD_HANDLE = @(x,varargin)zeros(size(x,1),1);  % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap;                      % Exact siolution
  TOL = .2;                                     % Stopping criterion
  THETA = .3;                                   % Refinement strategy
  NITER = 5;                                    % Plotting levels
  PLOTTING = 0;                                 % plotting flag
  zero_fun = @(x,varargin) zeros(size(x,1),1);  % left hand side for dual problem
  MAX_NDOFS = 10^4;
  fig = figure('Name','Output functional');
  
  limit = -1.8551;

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
      
  nDofs = 0;
  F=[];
  while(nDofs < MAX_NDOFS)
      
    nDofs = [nDofs size(Mesh.Coordinates,1)];  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    B = assemMat_LFE(Mesh,@MASS_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % compute finite element approximation of cutoff function phi
    
    [phi_interp,FreeDofs] = assemDir_LFE(Mesh,-1,@gd_phi);
    M = A+B;
    L = -M*phi_interp;
    phi_interp(FreeDofs) = M(FreeDofs,FreeDofs)\L(FreeDofs);
    
    
    % compute the solution of the dual problem
    
    L_d = A'*phi_interp;
    [gF,FreeDofs]=assemDir_LFE(Mesh,-1,zero_fun);
    L_d = L_d-A*gF;
    A_T = A';
    gF(FreeDofs) = A_T(FreeDofs,FreeDofs)\L_d(FreeDofs);
    
    Mesh = add_Patches(Mesh);
    Eta_K = ErrEst_GOAL(gF,U,F_HANDLE,Mesh,P7O6(),QuadRule_1D);
    
    Eta = sqrt(sum(Eta_K.^2));
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % compute output functional
      
    F = [F U'*A*phi_interp - integr(Mesh,phi_interp,P7O6)];
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    %Mesh = refine_LEB(Mesh,MarkedElem);
    Mesh = refine_REG(Mesh);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1; 
    
    % Print out information
    
    fprintf('Eta  :  %2.4e,  F(u-u_n)  :  %2.4e\n',Eta,abs(F(end)-limit));
    
    % update figure
    clf;
    subplot(2,1,1),plot(nDofs(2:end),F);
    set(gca,'XScale','log','YScale','lin');    
    subplot(2,1,2),plot(nDofs(2:end),abs(F-limit));
    set(gca,'XScale','log','YScale','log');
    drawnow;
  end
 
  
  plot(nDofs(2:end),F);
  set(gca,'XScale','log','YScale','lin');
  
  limit = F(end);
  
  