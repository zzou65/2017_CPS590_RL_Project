function Mloc = MASS_W1F2nd(Vertices,ElemInfo,MU_HANDLE,QuadRule,varargin)
% MASS_W1F element mass matrix with weight mu for edge elements in 2D
%
%   MLOC = MASS_W1F(VERTICES) computes the element mass matrix using 
%   Whitney 1-forms finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
% 
%   ElemInfo (not used)
%
%   MU_HANDLE handle to a functions expecting a matrix whose rows
%   represent position arguments. Return value must be a vector
%   (variable arguments will be passed to this function)
%
%   Example:
%
%   Mloc = MASS_W1F(Vertices,ElemInfo,MU_HANDLE,QuadRule);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  
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
  Fval = MU_HANDLE(x,ElemInfo,varargin{:});

  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
  b4 = baryc(:,2)*g3'+baryc(:,3)*g2';
  b5 = baryc(:,3)*g1'+baryc(:,1)*g3';
  b6 = baryc(:,1)*g2'+baryc(:,2)*g1';
    
  % Compute local mass matrix
  
  weights = QuadRule.w * det_BK;
  Mloc(1,1) = sum(weights.*Fval.*sum(b1.*b1,2));
  Mloc(2,2) = sum(weights.*Fval.*sum(b2.*b2,2));
  Mloc(3,3) = sum(weights.*Fval.*sum(b3.*b3,2));
  Mloc(4,4) = sum(weights.*Fval.*sum(b4.*b4,2));
  Mloc(5,5) = sum(weights.*Fval.*sum(b5.*b5,2));
  Mloc(6,6) = sum(weights.*Fval.*sum(b6.*b6,2));
  
  Mloc(1,2) = sum(weights.*Fval.*sum(b1.*b2,2)); Mloc(2,1) = Mloc(1,2);
  Mloc(1,3) = sum(weights.*Fval.*sum(b1.*b3,2)); Mloc(3,1) = Mloc(1,3);
  Mloc(1,4) = sum(weights.*Fval.*sum(b1.*b4,2)); Mloc(4,1) = Mloc(1,4);
  Mloc(1,5) = sum(weights.*Fval.*sum(b1.*b5,2)); Mloc(5,1) = Mloc(1,5);
  Mloc(1,6) = sum(weights.*Fval.*sum(b1.*b6,2)); Mloc(6,1) = Mloc(1,6);
  
  Mloc(2,3) = sum(weights.*Fval.*sum(b2.*b3,2)); Mloc(3,2) = Mloc(2,3);
  Mloc(2,4) = sum(weights.*Fval.*sum(b2.*b4,2)); Mloc(4,2) = Mloc(2,4);
  Mloc(2,5) = sum(weights.*Fval.*sum(b2.*b5,2)); Mloc(5,2) = Mloc(2,5);
  Mloc(2,6) = sum(weights.*Fval.*sum(b2.*b6,2)); Mloc(6,2) = Mloc(2,6);
  
  Mloc(3,4) = sum(weights.*Fval.*sum(b3.*b4,2)); Mloc(4,3) = Mloc(3,4);
  Mloc(3,5) = sum(weights.*Fval.*sum(b3.*b5,2)); Mloc(5,3) = Mloc(3,5);
  Mloc(3,6) = sum(weights.*Fval.*sum(b3.*b6,2)); Mloc(6,3) = Mloc(3,6);

  Mloc(4,5) = sum(weights.*Fval.*sum(b4.*b5,2)); Mloc(5,4) = Mloc(4,5);
  Mloc(4,6) = sum(weights.*Fval.*sum(b4.*b6,2)); Mloc(6,4) = Mloc(4,6);
  
  Mloc(5,6) = sum(weights.*Fval.*sum(b5.*b6,2)); Mloc(6,5) = Mloc(5,6);
  
  return