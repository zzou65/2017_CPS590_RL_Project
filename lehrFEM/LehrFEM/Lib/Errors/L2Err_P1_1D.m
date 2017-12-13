function err = L2Err_P1_1D(Coordinates,u,QuadRule,FHandle,varargin)
% L2ERR_P1_1D Discretization error in L2 norm for 1D linear finite elements.
%
%   ERR = L2ERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE) computes the value of
%   the discretization error between the exact solution given by the
%   function handle FHANDLE and the finite element solution U.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = L2ERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = L2Err_P1_1D(Coordinates,u,QuadRule,FHandle);
%
%   See also shap_P1_1D.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  
  nCoordinates = size(Coordinates,1);
  
  % Precompute shape function values at thhe quedrature points
  
  N = shap_P1_1D(QuadRule.x);
    
  % Compute discretization error

  err = 0;
  for i = 1:(nCoordinates-1)
       
    % Compute element mapping  
        
    h = abs(Coordinates(i+1)-Coordinates(i));
    x = Coordinates(i)+h*QuadRule.x;  
    
    % Evaluate solutions
      
    u_EX = FHandle(x,varargin{:});
    u_FE = u(i)*N(:,1)+u(i+1)*N(:,2);
      
    % Compute error on current element
      
    err = err+sum(QuadRule.w.*abs(u_EX-u_FE).^2)*h;
      
  end
  err = sqrt(err);
  
return