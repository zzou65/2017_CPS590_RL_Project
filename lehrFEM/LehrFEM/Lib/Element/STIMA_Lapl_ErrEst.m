function Diag = STIMA_Lapl_ErrEst(EdgeLoc,Vertices,ElemFlag,QuadRule,varargin)
% STIMA_LAPL_ERREST Diagonal entry of element stiffness matrix.
%
%   DIAG = STIMA_LAPL_ERREST(EDGELOC,VERTICES,ELEMFLAG,QUADRULE) computes
%   the diagonal entry of the element stiffness matrix for the Laplace
%   equation.
%
%   Example:
%
%   Diag = STIMA_Lapl_ErrEst(EdgeLoc,Vertices,ElemFlag,P7O6());

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);

  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
    
  % Compute values of gradient of test function

  grad_N = grad_shap_EP2(QuadRule.x);
  switch(EdgeLoc)
    case 1
      grad_N = grad_N(:,1:2)*transpose(inv_BK);      
    case 2
      grad_N = grad_N(:,3:4)*transpose(inv_BK);    
    case 3
      grad_N = grad_N(:,5:6)*transpose(inv_BK);    
  end
  
  % Compute diagonal entry of stiffness matrix on the current element
  
  Diag = sum(QuadRule.w.*sum(grad_N.*grad_N,2))*det_BK;
  
return
