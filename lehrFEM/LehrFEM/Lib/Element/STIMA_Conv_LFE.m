function Mloc =STIMA_Conv_LFE(Vertices, dummy ,VHANDLE, QuadRule, varargin)
% STIMA_Conv_LFE Element matrix for convection terms.
%
%   Mloc = STIMA_Conv_LFE(VERTICES, dummy, VHandle, QuadRule) computes the 
%   element matrix for convection using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   VHANDLE function handle for convection terms
%
%   Example:
%
%   Mloc = MASS_LFE(Vertices, d, Vhandle, P7O6());

%   Copyright 2005-2007 Patrick Meury, Holger Heumann, Alan Mitchell
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
 
% Calculate determinant of BK = twice the area of a triangle
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK=inv(BK);
  det_BK = abs(det(BK));

% Define x values

  nPoints = size(QuadRule.w,1); 
  x = QuadRule.x*BK + ones(nPoints,1)*bK;

% Compute velocity and shap functions at quadrature points
  
  v=VHANDLE(x,varargin{:});
  
  Lambda = shap_LFE(QuadRule.x);
  
  grad_lambda = grad_shap_LFE(QuadRule.x);
  
  grad_lambda(:,1:2)=grad_lambda(:,1:2)*transpose(inv_BK);
  grad_lambda(:,3:4)=grad_lambda(:,3:4)*transpose(inv_BK);
  grad_lambda(:,5:6)=grad_lambda(:,5:6)*transpose(inv_BK);
  
% Compute components of the local element matrix using linear
% Lagrangian finite elements
  
  Mloc(1,1) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,1:2).*v,2).*Lambda(:,1));
  Mloc(1,2) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,1:2).*v,2).*Lambda(:,2));
  Mloc(1,3) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,1:2).*v,2).*Lambda(:,3));
 
  Mloc(2,1) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,3:4).*v,2).*Lambda(:,1));
  Mloc(2,2) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,3:4).*v,2).*Lambda(:,2));
  Mloc(2,3) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,3:4).*v,2).*Lambda(:,3));
  
  Mloc(3,1) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,5:6).*v,2).*Lambda(:,1));
  Mloc(3,2) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,5:6).*v,2).*Lambda(:,2));
  Mloc(3,3) = det_BK*sum(QuadRule.w.*sum(grad_lambda(:,5:6).*v,2).*Lambda(:,3));
  
  Mloc=Mloc';
   
return