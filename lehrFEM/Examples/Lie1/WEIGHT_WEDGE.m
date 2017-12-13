function  [Mloc1 Mloc2 Mloc3] = WEIGHT_WEDGE(Vertices)
% MASS_W1F element mass matrices for edge elements in 2D on intersection 
% of triangle with cells
%
%   MLOC = MASS_W1F(VERTICES) computes the element mass matrices using 
%   Whitney 1-forms finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
% 
%
%   Example:
%
%   Mloc = WEIGHT_WEDGE(Vertices);

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  MU_HANDLE=@(x,varargin)1;

  QuadRule=P3O3();
  nPoints = size(QuadRule.w,1);
  % Compute basis functions on Triangle 
  
  % Compute element mapping  

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  m1 = (Vertices(2,:)+Vertices(3,:))/2;
  m2 = (Vertices(3,:)+Vertices(1,:))/2;
  m3 = (Vertices(1,:)+Vertices(2,:))/2;
  b =(Vertices(1,:)+Vertices(2,:)+Vertices(3,:))/3;
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  
  % 1. vertice 
  % a) Mass on first subtriangle <P1, b, m2>
   
  BK_h = [ b - P1 ; m2 - P1 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P1;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  % Compute local mass matrix

  weights = QuadRule.w * det_BK_h;
  Mloc1(1,1) = sum(weights.*sum(b1.*b1,2));
  Mloc1(2,2) = sum(weights.*sum(b2.*b2,2));
  Mloc1(3,3) = sum(weights.*sum(b3.*b3,2));
  Mloc1(1,2) = sum(weights.*sum(b1.*b2,2)); 
  Mloc1(1,3) = sum(weights.*sum(b1.*b3,2)); 
  Mloc1(2,3) = sum(weights.*sum(b2.*b3,2)); 
  
  % b) Mass on second subtriangle <P1, m3, b>
   
  BK_h = [ m3 - P1 ; b - P1 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P1;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  % Compute local mass matrix

  weights = QuadRule.w * det_BK_h;
  Mloc1(1,1) = Mloc1(1,1)+sum(weights.*sum(b1.*b1,2));
  Mloc1(2,2) = Mloc1(2,2)+sum(weights.*sum(b2.*b2,2));
  Mloc1(3,3) = Mloc1(3,3)+sum(weights.*sum(b3.*b3,2));
  Mloc1(1,2) = Mloc1(1,2)+sum(weights.*sum(b1.*b2,2)); Mloc1(2,1) = Mloc1(1,2);
  Mloc1(1,3) = Mloc1(1,3)+sum(weights.*sum(b1.*b3,2)); Mloc1(3,1) = Mloc1(1,3);
  Mloc1(2,3) = Mloc1(2,3)+sum(weights.*sum(b2.*b3,2)); Mloc1(3,2) = Mloc1(2,3);  
  
  % 2. vertex 
  % a) Mass on first subtriangle <P2, b, m3>
   
  BK_h = [ b - P2 ; m3 - P2 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P2;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  % Compute local mass matrix

  weights = QuadRule.w * det_BK_h;
  Mloc2(1,1) = sum(weights.*sum(b1.*b1,2));
  Mloc2(2,2) = sum(weights.*sum(b2.*b2,2));
  Mloc2(3,3) = sum(weights.*sum(b3.*b3,2));
  Mloc2(1,2) = sum(weights.*sum(b1.*b2,2)); 
  Mloc2(1,3) = sum(weights.*sum(b1.*b3,2)); 
  Mloc2(2,3) = sum(weights.*sum(b2.*b3,2)); 
  
  % b) Mass on first subtriangle <P2, m1, b>
   
  BK_h = [ m1 - P2 ; b - P2 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P2;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  % Compute local mass matrix
  weights = QuadRule.w * det_BK_h;
  Mloc2(1,1) = Mloc2(1,1)+sum(weights.*sum(b1.*b1,2));
  Mloc2(2,2) = Mloc2(2,2)+sum(weights.*sum(b2.*b2,2));
  Mloc2(3,3) = Mloc2(3,3)+sum(weights.*sum(b3.*b3,2));
  Mloc2(1,2) = Mloc2(1,2)+sum(weights.*sum(b1.*b2,2)); Mloc2(2,1) = Mloc2(1,2);
  Mloc2(1,3) = Mloc2(1,3)+sum(weights.*sum(b1.*b3,2)); Mloc2(3,1) = Mloc2(1,3);
  Mloc2(2,3) = Mloc2(2,3)+sum(weights.*sum(b2.*b3,2)); Mloc2(3,2) = Mloc2(2,3);
  
  % 3. vertice 
  % a) Mass on first subtriangle <P3, b, m1>
   
  BK_h = [ b - P3 ; m1 - P3 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P3;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  % Compute local mass matrix

  weights = QuadRule.w * det_BK_h;
  Mloc3(1,1) = sum(weights.*sum(b1.*b1,2));
  Mloc3(2,2) = sum(weights.*sum(b2.*b2,2));
  Mloc3(3,3) = sum(weights.*sum(b3.*b3,2));
  Mloc3(1,2) = sum(weights.*sum(b1.*b2,2)); 
  Mloc3(1,3) = sum(weights.*sum(b1.*b3,2)); 
  Mloc3(2,3) = sum(weights.*sum(b2.*b3,2)); 
  
  % b) Mass on second subtriangle <P3, m2, b>
   
  BK_h = [ m2 - P3 ; b - P3 ]; % transpose of transformation matrix
  det_BK_h = abs(det(BK_h));     % twice the area of the triagle
   
  x = QuadRule.x*BK_h + ones(nPoints,1)*P3;
  
  b1 = 1/det_BK*(x-ones(nPoints,1)*P1);
  b2 = 1/det_BK*(x-ones(nPoints,1)*P2);
  b3 = 1/det_BK*(x-ones(nPoints,1)*P3);
    
  
  % Compute local mass matrix
  weights = QuadRule.w * det_BK_h;
  Mloc3(1,1) = Mloc3(1,1)+sum(weights.*sum(b1.*b1,2));
  Mloc3(2,2) = Mloc3(2,2)+sum(weights.*sum(b2.*b2,2));
  Mloc3(3,3) = Mloc3(3,3)+sum(weights.*sum(b3.*b3,2));
  Mloc3(1,2) = Mloc3(1,2)+sum(weights.*sum(b1.*b2,2)); Mloc3(2,1) = Mloc3(1,2);
  Mloc3(1,3) = Mloc3(1,3)+sum(weights.*sum(b1.*b3,2)); Mloc3(3,1) = Mloc3(1,3);
  Mloc3(2,3) = Mloc3(2,3)+sum(weights.*sum(b2.*b3,2)); Mloc3(3,2) = Mloc3(2,3);
    
  
return