function Aloc = STIMA_Heat_LFE(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% STIMA_HEAT_LFE Element stiffness matrix for the stationary heat equation.
%
%   ALOC = STIMA_HEAT_QFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element stiffness matrix for the data given by function handle FHANDLE. 
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
%   FHANDLE is the function handle to the heat conductivity of the
%   stationary heat conduction equation.
%   
%   ALOC = STIMA_HEAT_LFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also
%   handles the additional variable length argument list FPARAM to the
%   function handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)1;   
%   QuadRule = P7O6_rule();
%   Aloc = STIMA_Heat_LFE([0 0; 1 0; 0 1],0,QuadRule,FHandle);
%
%   See also grad_shap_LFE.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPoints = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(3,3);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
  
  TK = det_BK*transpose(inv_BK)*inv_BK;
  
  x = QuadRule.x*BK+ones(nPoints,1)*bK;
  
  % Compute element stiffness matrix
  
  FVal = FHandle(x,ElemInfo,varargin{:});
  grad_N = grad_shap_LFE(QuadRule.x);
  
  Aloc(1,1) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,1:2)*TK),2));
  Aloc(1,2) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,3:4)*TK),2));
  Aloc(1,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,1:2)).*(grad_N(:,5:6)*TK),2));
  Aloc(2,2) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,3:4)*TK),2));
  Aloc(2,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,3:4)).*(grad_N(:,5:6)*TK),2));
  Aloc(3,3) = sum(QuadRule.w.*FVal.*sum((grad_N(:,5:6)).*(grad_N(:,5:6)*TK),2));
  
  % Update lower triangular part
  
  Aloc(2,1) = Aloc(1,2);
  Aloc(3,1) = Aloc(1,3);
  Aloc(3,2) = Aloc(2,3);
  
return