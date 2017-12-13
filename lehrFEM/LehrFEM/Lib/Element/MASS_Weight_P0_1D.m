function Mloc = MASS_Weight_P0_1D(Vertices,QuadRule,FHandle,varargin)
% MASS_WEIGHT_P0_1D Element mass matrix.
%
%   MLOC = MASS_WEIGHT_P0_1D(VERTICES) computes the element mass matrix
%   using constant finite elements.
%
%   VERTICES is 2-by-1 matrix specifying the vertices of the current
%   element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%    
%   FHANDLE is the function handle to the weight function.
%   
%   Example:
%
%   Mloc = MASS_Weight_P0_1D([0; 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping
  
  h = abs(Vertices(2)-Vertices(1));
  x = Vertices(1)+QuadRule.x*h;
  
  % Compute element mass matrix
  
  FVal = FHandle(x,varargin{:});  
  
  Mloc = sum(QuadRule.w.*FVal)*h;

return