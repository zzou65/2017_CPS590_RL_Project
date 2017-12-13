function Mloc = STIMA_GradContr2(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_  not working
%
%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));      % twice the area of the triagle
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points 
  nPoints = size(QuadRule.w,1);
  baryc= [QuadRule.x,1-sum(QuadRule.x,2)];

  % Quadrature points in actual element
  % stored as rows of a matrix
  x = QuadRule.x*BK + ones(nPoints,1)*P1;

  % Evaluate coefficient function at quadrature nodes
  Fval = V_HANDLE(1/3*sum(Vertices,1),ElemInfo,varargin{:});
  
  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
    
  % auxilary product: Gradient of linear shape functions times edges
  % element basis functions
  g1b1=b1*g1;
  g2b1=b1*g2;
  g3b1=b1*g3;
  
  g1b2=b2*g1;
  g2b2=b2*g2;
  g3b2=b2*g3;
  
  g1b3=b3*g1;
  g2b3=b3*g2;
  g3b3=b3*g3;
  
  % Compute local mass matrix
  weights = QuadRule.w * det_BK;
  Mloc(1,1) = Fval*g3*sum(weights.*g2b1)-Fval*g2*sum(weights.*g3b1);
  Mloc(1,2) = Fval*g1*sum(weights.*g3b1)-Fval*g3*sum(weights.*g1b1);
  Mloc(1,3) = Fval*g2*sum(weights.*g1b1)-Fval*g1*sum(weights.*g2b1);
  
  Mloc(2,1) = Fval*g3*sum(weights.*g2b2)-Fval*g2*sum(weights.*g3b2);
  Mloc(2,2) = Fval*g1*sum(weights.*g3b2)-Fval*g3*sum(weights.*g1b2);
  Mloc(2,3) = Fval*g2*sum(weights.*g1b2)-Fval*g1*sum(weights.*g2b2);
  
  Mloc(3,1) = Fval*g3*sum(weights.*g2b3)-Fval*g2*sum(weights.*g3b3);
  Mloc(3,2) = Fval*g1*sum(weights.*g3b3)-Fval*g3*sum(weights.*g1b3);
  Mloc(3,3) = Fval*g2*sum(weights.*g1b3)-Fval*g1*sum(weights.*g2b3);
  
  return