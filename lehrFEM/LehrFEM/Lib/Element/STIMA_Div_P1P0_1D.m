function Aloc = STIMA_Div_P1P0_1D(Vertices,QuadRule,FHandle,varargin)
% STIMA_DIV_P1P0_1D Element stiffness matrix.
%
%   ALOC = STIMA_DIV_P1P0_1D(VERTICES,QUADRULE,FHANDLE) computes the
%   element stiffness matrix for the data given by function handle FHANDLE. 
%
%   VERTICES is a 2-by-1 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%    
%   FHANDLE is the function handle to the weight function.
%   
%   ALOC = STIMA_HEAT_P1_1D(VERTICES,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function
%   handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)1;   
%   Aloc = STIMA_Div_P1P0_1D([0 0; 1 0; 0 1],gauleg(0,1,5),FHandle);
%
%   See also grad_shap_P1_1D.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Preallocate memory
  
  Aloc = zeros(2,1);
  
  % Compute element mapping
  
  h = abs(Vertices(2)-Vertices(1));
  x = Vertices(1)+QuadRule.x*h;
  
  % Compute element stiffness matrix
  
  FVal = FHandle(x,varargin{:});
  grad_N = grad_shap_P1_1D(QuadRule.x);
  
  Aloc(1) = sum(QuadRule.w.*FVal.*grad_N(:,1));
  Aloc(2) = sum(QuadRule.w.*FVal.*grad_N(:,2));

return