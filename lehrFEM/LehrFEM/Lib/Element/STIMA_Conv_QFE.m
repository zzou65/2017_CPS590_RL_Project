function Mloc =STIMA_Conv_QFE(Vertices, dummy ,VHANDLE, QuadRule, Varargin)
% STIMA_Conv_QFE Element matrix for convection terms.
%
%   Mloc = STIMA_Conv_QFE(VERTICES, dummy, VHandle, QuadRule) computes the 
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

%   Copyright 2005-2007 Patrick Meury, Holger Heumann
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
  
  v=VHANDLE(x);
  
  Lambda = shap_QFE(QuadRule.x);
  
  grad_lambda = grad_shap_QFE(QuadRule.x);
  
  grad_lambda(:,1:2)=grad_lambda(:,1:2)*transpose(inv_BK);
  grad_lambda(:,3:4)=grad_lambda(:,3:4)*transpose(inv_BK);
  grad_lambda(:,5:6)=grad_lambda(:,5:6)*transpose(inv_BK);
  grad_lambda(:,7:8)=grad_lambda(:,7:8)*transpose(inv_BK);
  grad_lambda(:,9:10)=grad_lambda(:,9:10)*transpose(inv_BK);
  grad_lambda(:,11:12)=grad_lambda(:,11:12)*transpose(inv_BK);
  
% Compute components of the local element matrix using bilinear
% Lagrangian finite elements
  
  vgl=sum(grad_lambda(:,1:2).*v,2);
  Mloc(1,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(1,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(1,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(1,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(1,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(1,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  vgl=sum(grad_lambda(:,3:4).*v,2);
  Mloc(2,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(2,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(2,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(2,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(2,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(2,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  vgl=sum(grad_lambda(:,5:6).*v,2);
  Mloc(3,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(3,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(3,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(3,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(3,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(3,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  vgl=sum(grad_lambda(:,7:8).*v,2);
  Mloc(4,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(4,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(4,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(4,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(4,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(4,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  vgl=sum(grad_lambda(:,9:10).*v,2);
  Mloc(5,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(5,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(5,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(5,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(5,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(5,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  vgl=sum(grad_lambda(:,11:12).*v,2);
  Mloc(6,1) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,1));
  Mloc(6,2) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,2));
  Mloc(6,3) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,3));
  Mloc(6,4) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,4));
  Mloc(6,5) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,5));
  Mloc(6,6) = det_BK*sum(QuadRule.w.*vgl.*Lambda(:,6));
  
  Mloc=Mloc';
   
return