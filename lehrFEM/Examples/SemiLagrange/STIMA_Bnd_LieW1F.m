function Aloc = STIMA_Bnd_LieW1F(Edge,Normal,BdFlag,Data,V_Handle, QuadRule,varargin)
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
  
  Aloc = zeros(3,3);
  nPoints = size(QuadRule.w,1);
 
  dS = norm(Edge(2,:)-Edge(1,:));
  x = QuadRule.x*(Edge(2,:)-Edge(1,:)) + ones(nPoints,1)*Edge(1,:);
  Fval= V_Handle(x);

  % Compute on the left and right hand side element
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:)-bK; ...
            Data.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK=abs(det(BK));
  
  % Compute constant gradients of barycentric coordinate functions
  g1 = [Data.Vertices(2,2)-Data.Vertices(3,2);Data.Vertices(3,1)-Data.Vertices(2,1)]/det_BK;
  g2 = [Data.Vertices(3,2)-Data.Vertices(1,2);Data.Vertices(1,1)-Data.Vertices(3,1)]/det_BK;
  g3 = [Data.Vertices(1,2)-Data.Vertices(2,2);Data.Vertices(2,1)-Data.Vertices(1,1)]/det_BK;
  
  % Get barycentric coordinates of quadrature points
  if Data.Match
      switch(Data.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) 1-QuadRule.x QuadRule.x];
          case 2
              baryc = [QuadRule.x zeros(nPoints,1) 1-QuadRule.x];
          case 3
              baryc = [1-QuadRule.x QuadRule.x zeros(nPoints,1)];
      end
  else
      switch(Data.EdgeLoc)
          case 1
              baryc = [zeros(nPoints,1) QuadRule.x 1-QuadRule.x];
          case 2
              baryc = [1-QuadRule.x zeros(nPoints,1) QuadRule.x];
          case 3
              baryc = [QuadRule.x 1-QuadRule.x zeros(nPoints,1)];
      end
  end
  N1 = baryc(:,2)*g3'-baryc(:,3)*g2';
  N2 = baryc(:,3)*g1'-baryc(:,1)*g3';
  N3 = baryc(:,1)*g2'-baryc(:,2)*g1';
  
  V_N1=sum(N1.*Fval,2);
  V_N2=sum(N2.*Fval,2);
  V_N3=sum(N3.*Fval,2);

  N_N1=N1*Normal';
  N_N2=N2*Normal';
  N_N3=N3*Normal';

  % Compute entries of element matrix
  
  gamma=dS;
  Aloc(1,1) = gamma*sum(QuadRule.w.*N_N1.*V_N1);
  Aloc(1,2) = gamma*sum(QuadRule.w.*N_N1.*V_N2);
  Aloc(1,3) = gamma*sum(QuadRule.w.*N_N1.*V_N3);
 
  Aloc(2,1) = gamma*sum(QuadRule.w.*N_N2.*V_N1);
  Aloc(2,2) = gamma*sum(QuadRule.w.*N_N2.*V_N2);
  Aloc(2,3) = gamma*sum(QuadRule.w.*N_N2.*V_N3);
  
  Aloc(3,1) = gamma*sum(QuadRule.w.*N_N3.*V_N1);
  Aloc(3,2) = gamma*sum(QuadRule.w.*N_N3.*V_N2);
  Aloc(3,3) = gamma*sum(QuadRule.w.*N_N3.*V_N3);
  
return