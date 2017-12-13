function A = STIMA_Div_y_P1P0(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% STIMA_DIV_Y_P1P0 Element stiffness matrix.
%
%   A = STIMA_DIV_Y_P1P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element stiffness matrix using linear/constant Lagrangian finite
%   elements.
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
%   A = STIMA_DIV_Y_P1P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   A = STIMA_Div_y_P1P0([0 0; 1 0; 0 1],0,P7O6(),FHandle);
%
%   See also grad_shap_LFE.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);

  % Preallocate memory
  
  A = zeros(3,1);

  % Compute element map  

  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK = abs(det(BK));
  
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Evaluate function
  
  FVal = FHandle(x,ElemInfo,varargin{:});
  
  % Compute element shape functions

  grad_N = grad_shap_LFE(QuadRule.x);
  grad_N(:,1:2) = grad_N(:,1:2)*inv_BK_t;
  grad_N(:,3:4) = grad_N(:,3:4)*inv_BK_t;
  grad_N(:,5:6) = grad_N(:,5:6)*inv_BK_t;
  
  % Compute entries of stiffness matrix
  
  A(1) = sum(QuadRule.w.*FVal.*grad_N(:,2))*det_BK; 
  A(2) = sum(QuadRule.w.*FVal.*grad_N(:,4))*det_BK;
  A(3) = sum(QuadRule.w.*FVal.*grad_N(:,6))*det_BK;
  
return