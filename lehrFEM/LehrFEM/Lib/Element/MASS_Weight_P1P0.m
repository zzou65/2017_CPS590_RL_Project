function M = MASS_Weight_P1P0(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% MASS_WEIGHT_P1P0 Element mass matrix.
%
%   M = MASS_WEIGHT_P1P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element mass matrix using linear/constant Lagrangian finite elements.
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
%   FHANDLE denotes the function handle to weight function.
%
%   M = MASS_WEIGHT_P1P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   M = MASS_Weight_P1P0([0 0; 1 0; 0 1],0,P7O6(),FHandle);
%
%   See also shap_LFE.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);

  % Preallocate memory
  
  M = zeros(3,1);
  
  % Evaluate shape functions
  
  N = shap_LFE(QuadRule.x);

  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [bK-Vertices(2,:); bK-Vertices(3,:)];
  det_BK = abs(det(BK));
  
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Evaluate function
  
  Fval = FHandle(x,ElemInfo,varargin);
  
  % Compute the matrix
  
  for i = 1:3
    M(i) = sum(QuadRule.w.*Fval.*N(:,i))*det_BK; 
  end
  
return