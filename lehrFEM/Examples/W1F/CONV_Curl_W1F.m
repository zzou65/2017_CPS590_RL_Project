function Aloc = CONV_Curl_W1F(Vertices,ElemInfo,W_HANDLE,QuadRule,varargin)
% CONV_Curl_W1F convection term with velocity field W_HANDLE

%   Copyright 2009-2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  BBK = inv(BK);
  
  % quadrature points on the actual element and on reference element
  nPoints = size(QuadRule.w,1);
  xi = QuadRule.x;
  x = QuadRule.x*BK + ones(nPoints,1)*P1;
  
  % compute the velocity field at the quadrature points
  w_val = W_HANDLE(x);
  cross_w(:,1) = w_val(:,2); 
  cross_w(:,2) = -w_val(:,1); 
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
    
  % Get barycentric coordinates of quadrature points 
  baryc= [QuadRule.x,1-sum(QuadRule.x,2)];
  
  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';

  % use the quadrature rule to integrate over the reference triangle
  weights = QuadRule.w;
  Aloc = zeros(3,3);
  Aloc(1,1)=sum(weights.*sum(2*cross_w.*b1,2));
  Aloc(2,1)=Aloc(1,1);
  Aloc(3,1)=Aloc(1,1);
  Aloc(1,2)=sum(weights.*sum(cross_w.*b2,2))*2;
  Aloc(2,2)=Aloc(1,2);
  Aloc(3,2)=Aloc(1,2);
  Aloc(1,3)=sum(weights.*sum(cross_w.*b3,2))*2;
  Aloc(2,3)=Aloc(1,3);
  Aloc(3,3)=Aloc(1,3);
  
  Aloc = Aloc';
  return
  
  