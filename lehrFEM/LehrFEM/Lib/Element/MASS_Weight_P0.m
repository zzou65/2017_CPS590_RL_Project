function M = MASS_Weight_P0(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% MASS_HEAT_P0 Element mass matrix.
%
%   M = MASS_WEIGHT_P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element mass matrix using constant Lagrangian finite elements.
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   FHANDLE denotes the function handle to the weight functions.
%
%   M = MASS_WEIGHT_P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   M = MASS_Weight_P0([0 0; 1 0; 0 1],0,P7O6(),FHandle);
%
%   See also MASS_BFE.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);

  % Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK ; Vertices(3,:)-bK];
  det_BK = abs(det(BK));
  
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Evaluate function
  
  Fval = FHandle(x,ElemInfo,varargin{:});
  
  % Compute the matrix
  
  M = sum(QuadRule.w.*Fval)*det_BK; 
  
return