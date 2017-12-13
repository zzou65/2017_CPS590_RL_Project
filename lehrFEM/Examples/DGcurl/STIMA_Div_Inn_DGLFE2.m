function Aloc = STIMA_Div_Inn_DGLFE2(Edge,Normal,BdFlag,LData,RData,s,varargin)
% STIMA_DIV_INN_DGLFE2 Element stiffness matrix for interior terms.
%
%   ALOC = DIV_INN_DGLFE2(EDGE,NORMAL,BDFLAG,LDATA,RDATA,S) computes the
%   entries of the element stiffness matrix for the interior terms for vector valued function.
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
%   The integer S can specifies wheter the diffusive fluxes are discretized
%   in a symmetric or anti-symmetric way:
%    +1 Antisymmetric discretization of diffusive fluxes
%    -1 Symmetric discretization of diffusive fluxes

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(12,12);
 
  % Compute normal derivatives on the left and right hand side element
  
  dS = norm(Edge(2,:)-Edge(1,:));
  
  
  LbK = LData.Vertices(1,:);
  LBK = [LData.Vertices(2,:)-LbK; ...
        LData.Vertices(3,:)-LbK];
  Larea = abs(det(LBK));
  
  INL = zeros(1,6);
  dNL = zeros(1,6);
  
  dNL = [ LData.Vertices(2,2) - LData.Vertices(3,2) ...
        LData.Vertices(3,1) - LData.Vertices(2,1) ...
        LData.Vertices(3,2) - LData.Vertices(1,2) ...
        LData.Vertices(1,1) - LData.Vertices(3,1) ...
        LData.Vertices(1,2) - LData.Vertices(2,2) ...
        LData.Vertices(2,1) - LData.Vertices(1,1) ]/(Larea);
 
  switch(LData.EdgeLoc)
    case 1
      INL(1) = 0;
      INL(2) = 0;
      INL(3) = Normal(1)*dS/2;
      INL(4) = Normal(2)*dS/2;
      INL(5) = Normal(1)*dS/2;
      INL(6) = Normal(2)*dS/2;
    case 2
      INL(1) = Normal(1)*dS/2;
      INL(2) = Normal(2)*dS/2;
      INL(3) = 0;
      INL(4) = 0;
      INL(5) = Normal(1)*dS/2;
      INL(6) = Normal(2)*dS/2;
    case 3
      INL(1) = Normal(1)*dS/2;
      INL(2) = Normal(2)*dS/2;
      INL(3) = Normal(1)*dS/2;
      INL(4) = Normal(2)*dS/2;
      INL(5) = 0;
      INL(6) = 0;
  end
  
  RbK = RData.Vertices(1,:);
  RBK = [RData.Vertices(2,:)-RbK; ...
        RData.Vertices(3,:)-RbK];
  Rarea = abs(det(RBK));
  
  INR = zeros(1,6);
  dNR = zeros(1,6);
  
  dNR = [ RData.Vertices(2,2) - RData.Vertices(3,2) ...
        RData.Vertices(3,1) - RData.Vertices(2,1) ...
        RData.Vertices(3,2) - RData.Vertices(1,2) ...
        RData.Vertices(1,1) - RData.Vertices(3,1) ...
        RData.Vertices(1,2) - RData.Vertices(2,2) ...
        RData.Vertices(2,1) - RData.Vertices(1,1) ]/(Rarea); 

  switch(RData.EdgeLoc)
    case 1
      INR(1) = 0;
      INR(2) = 0;
      INR(3) = Normal(1)*dS/2;
      INR(4) = Normal(2)*dS/2;
      INR(5) = Normal(1)*dS/2;
      INR(6) = Normal(2)*dS/2;
    case 2
      INR(1) = Normal(1)*dS/2;
      INR(2) = Normal(2)*dS/2;
      INR(3) = 0;
      INR(4) = 0;
      INR(5) = Normal(1)*dS/2;
      INR(6) = Normal(2)*dS/2;
    case 3
      INR(1) = Normal(1)*dS/2;
      INR(2) = Normal(2)*dS/2;
      INR(3) = Normal(1)*dS/2;
      INR(4) = Normal(2)*dS/2;
      INR(5) = 0;
      INR(6) = 0;
  end
    
  % Compute entries of element penalty matrix
  
  if(LData.Element < RData.Element)
    gamma = 1;
  else
    gamma = -1;
  end
  
  for i = 1:6
      for j = 1:6
        Aloc(i,j) = gamma*(dNL(i)*INL(j)+s*INL(i)*dNL(j))/2;
        Aloc(i,6+j) = gamma*(-dNL(i)*INR(j)+s*INL(i)*dNR(j))/2;
        Aloc(6+i,j) = gamma*(dNR(i)*INL(j)-s*INR(i)*dNL(j))/2;
        Aloc(6+i,6+j) = gamma*(-dNR(i)*INR(j)-s*INR(i)*dNR(j))/2;
      end
  end
 
return