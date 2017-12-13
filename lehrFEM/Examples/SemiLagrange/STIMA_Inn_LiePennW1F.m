function Aloc = STIMA_Inn_LieW1F(Edge,Normal,BdFlag,LData,RData,V_Handle, QuadRule,ce, varargin)
%   need to be tested  
% 
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   The integer BDFLAG denotes the boundary flag of the current edge. Note
%   that for interior edges only values larger than are allowed.
%
%   The structs LDATA and RDATA conatin the left and right hand side
%   element data:
%    ELEMENT  Integer specifying the neighbouring element.
%    ELEMFLAG Integer specifying the element flag of the neighbouring
%             element or zero.
%    VERTICES 3-by-2 matrix specifying the vertices of the neighbouring
%             element.
%    EDGELOC  Integer specifying the local edge number on the neighbouring
%             element.
%    MATCH    Integer specifying the relative orientation of the edge with
%             respect to the orientation of the neighbouring element.
%
%   See also grad_shap_DGCR.

%   Copyright 2006-2009 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(6,6);
  nPoints = size(QuadRule.w,1);
 
  dS = norm(Edge(2,:)-Edge(1,:));
  x = QuadRule.x*(Edge(2,:)-Edge(1,:)) + ones(nPoints,1)*Edge(1,:);
  FVal= V_Handle(x);

  % Compute on the left and right hand side element
   
  bK = LData.Vertices(1,:);
  BK = [LData.Vertices(2,:)-bK; ...
            LData.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK=abs(det(BK));
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [LData.Vertices(2,2)-LData.Vertices(3,2);LData.Vertices(3,1)-LData.Vertices(2,1)]/det_BK;
  g2 = [LData.Vertices(3,2)-LData.Vertices(1,2);LData.Vertices(1,1)-LData.Vertices(3,1)]/det_BK;
  g3 = [LData.Vertices(1,2)-LData.Vertices(2,2);LData.Vertices(2,1)-LData.Vertices(1,1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points
  if LData.Match==1
      switch(LData.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) 1-QuadRule.x QuadRule.x];
          case 2
              baryc = [QuadRule.x zeros(nPoints,1) 1-QuadRule.x];
          case 3
              baryc = [1-QuadRule.x QuadRule.x zeros(nPoints,1)];
      end
  else
      switch(LData.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) QuadRule.x 1-QuadRule.x];
          case 2
              baryc = [1-QuadRule.x zeros(nPoints,1) QuadRule.x];
          case 3
              baryc = [QuadRule.x 1-QuadRule.x zeros(nPoints,1)];
      end
  end
  NL1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  NL2= baryc(:,3)*g1'-baryc(:,1)*g3';
  NL3 = baryc(:,1)*g2'-baryc(:,2)*g1';
  
  N_NL1=NL1*Normal';
  N_NL2=NL2*Normal';
  N_NL3=NL3*Normal';
  
  bK = RData.Vertices(1,:);
  BK = [RData.Vertices(2,:)-bK; ...
            RData.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK=abs(det(BK));
 
  % Compute constant gradients of barycentric coordinate functions
  g1 = [RData.Vertices(2,2)-RData.Vertices(3,2);RData.Vertices(3,1)-RData.Vertices(2,1)]/det_BK;
  g2 = [RData.Vertices(3,2)-RData.Vertices(1,2);RData.Vertices(1,1)-RData.Vertices(3,1)]/det_BK;
  g3 = [RData.Vertices(1,2)-RData.Vertices(2,2);RData.Vertices(2,1)-RData.Vertices(1,1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points
  if RData.Match==1
      switch(RData.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) 1-QuadRule.x QuadRule.x];
          case 2
              baryc = [QuadRule.x zeros(nPoints,1) 1-QuadRule.x];
          case 3
              baryc = [1-QuadRule.x QuadRule.x  zeros(nPoints,1)];
      end
  else
      switch(RData.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) QuadRule.x 1-QuadRule.x];
          case 2
              baryc = [1-QuadRule.x zeros(nPoints,1) QuadRule.x];
          case 3
              baryc = [QuadRule.x 1-QuadRule.x zeros(nPoints,1)];
      end
  end

  NR1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  NR2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  NR3 = baryc(:,1)*g2'-baryc(:,2)*g1';
 
  N_NR1=NR1*Normal';
  N_NR2=NR2*Normal';
  N_NR3=NR3*Normal';
  
  % Compute entries of element matrix
  
%   if(LData.Element < RData.Element)
%     gamma = -abs(sum(FVal*Normal'))/nPoints;
%   else
%     gamma = -abs(sum(FVal*Normal'))/nPoints;
%   end

  gamma = abs(sum(FVal*Normal'))/nPoints;

  Aloc(1,1) = gamma*sum(QuadRule.w.*N_NL1.*N_NL1);
  Aloc(1,2) = gamma*sum(QuadRule.w.*N_NL1.*N_NL2);
  Aloc(1,3) = gamma*sum(QuadRule.w.*N_NL1.*N_NL3);
  
  Aloc(1,4) = - gamma*sum(QuadRule.w.*N_NL1.*N_NR1);
  Aloc(1,5) = - gamma*sum(QuadRule.w.*N_NL1.*N_NR2);
  Aloc(1,6) = - gamma*sum(QuadRule.w.*N_NL1.*N_NR3);

  Aloc(2,1) = gamma*sum(QuadRule.w.*N_NL2.*N_NL1);
  Aloc(2,2) = gamma*sum(QuadRule.w.*N_NL2.*N_NL2);
  Aloc(2,3) = gamma*sum(QuadRule.w.*N_NL2.*N_NL3);
  
  Aloc(2,4) = - gamma*sum(QuadRule.w.*N_NL2.*N_NR1);
  Aloc(2,5) = - gamma*sum(QuadRule.w.*N_NL2.*N_NR2);
  Aloc(2,6) = - gamma*sum(QuadRule.w.*N_NL2.*N_NR3);

  Aloc(3,1) = gamma*sum(QuadRule.w.*N_NL3.*N_NL1);
  Aloc(3,2) = gamma*sum(QuadRule.w.*N_NL3.*N_NL2);
  Aloc(3,3) = gamma*sum(QuadRule.w.*N_NL3.*N_NL3);
  
  Aloc(3,4) = - gamma*sum(QuadRule.w.*N_NL3.*N_NR1);
  Aloc(3,5) = - gamma*sum(QuadRule.w.*N_NL3.*N_NR2);
  Aloc(3,6) = - gamma*sum(QuadRule.w.*N_NL3.*N_NR3);

  Aloc(4,1) = - gamma*sum(QuadRule.w.*N_NR1.*N_NL1);
  Aloc(4,2) = - gamma*sum(QuadRule.w.*N_NR1.*N_NL2);
  Aloc(4,3) = - gamma*sum(QuadRule.w.*N_NR1.*N_NL3);
  
  Aloc(4,4) = gamma*sum(QuadRule.w.*N_NR1.*N_NR1);
  Aloc(4,5) = gamma*sum(QuadRule.w.*N_NR1.*N_NR2);
  Aloc(4,6) = gamma*sum(QuadRule.w.*N_NR1.*N_NR3);

  Aloc(5,1) = - gamma*sum(QuadRule.w.*N_NR2.*N_NL1);
  Aloc(5,2) = - gamma*sum(QuadRule.w.*N_NR2.*N_NL2);
  Aloc(5,3) = - gamma*sum(QuadRule.w.*N_NR2.*N_NL3);
  
  Aloc(5,4) = gamma*sum(QuadRule.w.*N_NR2.*N_NR1);
  Aloc(5,5) = gamma*sum(QuadRule.w.*N_NR2.*N_NR2);
  Aloc(5,6) = gamma*sum(QuadRule.w.*N_NR2.*N_NR3);

  Aloc(6,1) = - gamma*sum(QuadRule.w.*N_NR3.*N_NL1);
  Aloc(6,2) = - gamma*sum(QuadRule.w.*N_NR3.*N_NL2);
  Aloc(6,3) = - gamma*sum(QuadRule.w.*N_NR3.*N_NL3);
  
  Aloc(6,4) = gamma*sum(QuadRule.w.*N_NR3.*N_NR1);
  Aloc(6,5) = gamma*sum(QuadRule.w.*N_NR3.*N_NR2);
  Aloc(6,6) = gamma*sum(QuadRule.w.*N_NR3.*N_NR3);
  
  Aloc=ce*dS*Aloc;
  
return