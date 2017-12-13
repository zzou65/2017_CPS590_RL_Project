function err = H1SErr_P1_1D(Coordinates,u,QuadRule,FHandle,varargin)
% H1SERR_P1_1D Discretization error in H1 semi-norm for 1D linear finite
%            elements.
%
%   ERR = HISERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE) computes the value
%   of the discretization error between the exact solution given by the
%   function handle FHANDLE and the finite element solution U.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   ERR = HISERR_P1_1D(COORDINATES,U,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = HISErr_P1_1D(Coordinates,u,QuadRule,FHandle);
%
%   See also grad_shap_P1_1D.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,1);
    
  % Precompute gradients of shape functions
  
  grad_N = grad_shap_P1_1D(QuadRule.x);
  
  % Compute discretization error 
   
  err = 0;
  for i= 1:(nCoordinates-1)
          
    % Compute element mapping
      
    h = abs(Coordinates(i+1)-Coordinates(i));
    x = Coordinates(i)+h*QuadRule.x;
    
    % Evaluate solutions
      
    grad_u_EX = FHandle(x,varargin{:});
    grad_u_FE = (u(i)*grad_N(:,1)+u(i+1)*grad_N(:,2))/h;
      
    % Compute error on the current element
     
    err = err+sum(QuadRule.w.*abs(grad_u_FE-grad_u_EX).^2)*h;
      
  end
  err = sqrt(err);
    
return