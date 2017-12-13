function Mloc = STIMA_ContrRot_UPQuad(Vertices,ElemInfo,V_HANDLE,varargin)
% STIMA_ContrRot computes element contribution of -v x curl u term in
% standard FEM-variational setting
%
%   MLOC = STIMA_ContrRot(VERTICES ...) computes element contribution of v x curl u term matrix using 
%   Whitney 1-forms finite elements.
%
%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathem8atics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  QuadRule=P302();
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [P2(2)-P3(2);P3(1)-P2(1)]/det_BK;
  g2 = [P3(2)-P1(2);P1(1)-P3(1)]/det_BK;
  g3 = [P1(2)-P2(2);P2(1)-P1(1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points 
  nPoints = size(QuadRule.w,1);
  baryc= [1-sum(QuadRule.x,2),QuadRule.x];

  % Quadrature points in actual element
  % stored as rows of a matrix
  x = QuadRule.x*BK + ones(nPoints,1)*P1;

  % Evaluate coefficient function at quadrature nodes
  Fval = V_HANDLE(x,ElemInfo);
  
  EdgeContr=varargin{:};
  
  %Rotation
  Fval =[Fval(:,2) -Fval(:,1)];

  % Evaluate basis functions at quadrature points
  % the rows of b(i) store the value of the the i-th
  % basis function at the quadrature points
  b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
  
  % upwind curls
  b1=EdgeContr(1,:);
  b2=EdgeContr(2,:);
  b3=EdgeContr(3,:);
    
  % Compute local mass matrix
  weights = QuadRule.w * det_BK;
  Mloc(1,1) = [sum(weights.*sum(Fval.*b1,2).db1)];
  Mloc(1,2) = [sum(weights.*sum(Fval.*b1,2).db2)];
  Mloc(1,3) = [sum(weights.*sum(Fval.*b1,2).db3)];
  
  Mloc(2,1) = [sum(weights.*sum(Fval.*b2,2).db1)];
  Mloc(2,2) = [sum(weights.*sum(Fval.*b2,2).db2)];
  Mloc(2,3) = [sum(weights.*sum(Fval.*b2,2).db3)];
  
  Mloc(3,1) = [sum(weights.*sum(Fval.*b2,2).db1)];
  Mloc(3,2) = [sum(weights.*sum(Fval.*b2,2).db2)];
  Mloc(3,3) = [sum(weights.*sum(Fval.*b2,2).db3)];
  
  Mloc=-Mloc;
  return