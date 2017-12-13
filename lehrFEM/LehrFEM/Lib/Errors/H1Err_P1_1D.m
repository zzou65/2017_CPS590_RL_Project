function err = H1Err_P1_1D(Coordinates,u,QuadRule,FHandle,varargin)
% H1ERR_P1_1D Discretization error in H1 norm for 1D linear finite elements.
%
%   ERR = H1ERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE) computes the value of
%   the discretization error between the exact solution given by the
%   function handle FHANDLE and the finite element solution U.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = H1ERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = H1Err_P1_1D(Mesh,u,QuadRule,FHandle);
%
%   See also shap_P1_1D, grad_shap_P1_1D.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  % Inititialize constants
  
  nCoordinates = size(Coordinates,1);
  
  % Precompute gradients and values of shape functions
  
  N = shap_P1_1D(QuadRule.x);
  grad_N = grad_shap_P1_1D(QuadRule.x);
    
  % Compute discretization error
  
  err = 0;  
  for i= 1:(nCoordinates-1)
        
    % Compute element mapping
      
    h = abs(Coordinates(i+1)-Coordinates(i));
    x = Coordinates(i)+h*QuadRule.x;
              
    % Evaluate solutions
      
    [u_EX,grad_u_EX] = FHandle(x,varargin{:});
    u_FE = u(i)*N(:,1)+u(i+1)*N(:,2);
    grad_u_FE = (u(i)*grad_N(:,1)+u(i+1)*grad_N(:,2))/h;
      
    % Compute error on the current element
      
    err = err+sum(QuadRule.w.*(abs(u_EX-u_FE).^2+abs(grad_u_FE-grad_u_EX).^2))*h;
      
  end 
  err = sqrt(err);
    
return