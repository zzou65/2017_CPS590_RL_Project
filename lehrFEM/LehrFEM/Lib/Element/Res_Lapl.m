function Res = Res_Lapl(U,EdgeLoc,Vertices,ElemFlag,QuadRule,FHandle,varargin)
% RES_LAPL Left and right hand-side residual for the current edge
%
%   RES = RES_LAPL(U,EDGELOC,VERTICES,ELEMFLAG,QUADRULE,FHANDLE) computes 
%   the residual for the Laplace equation on the current element.
%
%   Example:
%
%   Res = RHandle(U,EdgeLoc,Vertices,1,P7O6(),@f_LShap);

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
  
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Compute load function values

  FVal = FHandle(x,ElemFlag,varargin{:});

  % Compute values of gradient of numerical solution
  
  grad_N = grad_shap_LFE(QuadRule.x);
  grad_U = (U(1)*grad_N(:,1:2) + ...
            U(2)*grad_N(:,3:4) + ...
            U(3)*grad_N(:,5:6))*transpose(inv_BK);
  
  % Compute gradient of test function
        
  N = shap_EP2(QuadRule.x);
  grad_N = grad_shap_EP2(QuadRule.x);
  switch(EdgeLoc)
    case 1
      N = N(:,1);
      grad_N = grad_N(:,1:2)*transpose(inv_BK);
    case 2
      N = N(:,2);
      grad_N = grad_N(:,3:4)*transpose(inv_BK);    
    case 3
      N = N(:,3);
      grad_N = grad_N(:,5:6)*transpose(inv_BK);
  end
      
  % Compute residual on the current element
  
  Res = sum(QuadRule.w.*(FVal.*N - sum(grad_U.*grad_N,2)))*det_BK;
  
return