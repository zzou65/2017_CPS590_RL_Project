function Mloc = MASS_Weight_P1_1D(Vertices,QuadRule,FHandle,varargin)
% MASS_WEIGHT_P1_1D Element mass matrix.
%
%   MLOC = MASS_WEIGHT_P1_1D(VERTICES) computes the element mass matrix
%   using linear finite elements.
%
%   VERTICES is 2-by-1 matrix specifying the vertices of the current
%   element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%    
%   FHANDLE is the function handle to the weight functionn.
%   
%   Example:
%
%   Mloc = MASS_Weight_P1_1D([0; 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Preallocate memory
  
  Mloc = zeros(2,2);

  % Compute element mapping
  
  h = abs(Vertices(2)-Vertices(1));
  x = Vertices(1)+QuadRule.x*h;
  
  % Compute element mass matrix
  
  FVal = FHandle(x,varargin{:});  
  N = shap_P1_1D(QuadRule.x);
  
  Mloc(1,1) = sum(QuadRule.w.*FVal.*N(:,1).*N(:,1))*h;
  Mloc(1,2) = sum(QuadRule.w.*FVal.*N(:,1).*N(:,2))*h;
  Mloc(2,2) = sum(QuadRule.w.*FVal.*N(:,2).*N(:,2))*h;
  
  % Set lower triangular part
  
  Mloc(2,1) = Mloc(1,2);

return