function Aloc = STIMA_Grad_Inn_DGLFE2(Edge,Normal,BdFlag,LData,RData,s,varargin)
% STIMA_GRAD_NN_DGLFE2 Element stiffness matrix for interior terms.
%
%   ALOC = STIMA_GRAD_NN_DGLFE2(EDGE,NORMAL,BDFLAG,LDATA,RDATA,S) computes the
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
  
  grad_shap = grad_shap_DGLFE([0 0]);
  
  bK = LData.Vertices(1,:);
  BK = [LData.Vertices(2,:)-bK; ...
        LData.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  
  INL = zeros(1,3);
  dNL = zeros(1,3);
  switch(LData.EdgeLoc)
    case 1
      INL(1) = 0;
      INL(2) = dS/2;
      INL(3) = dS/2;
    case 2
      INL(1) = dS/2;
      INL(2) = 0;
      INL(3) = dS/2;
    case 3
      INL(1) = dS/2;
      INL(2) = dS/2;
      INL(3) = 0;
  end
  dNL(1) = sum((grad_shap(1:2)*inv_BK_t).*Normal,2);
  dNL(2) = sum((grad_shap(3:4)*inv_BK_t).*Normal,2);
  dNL(3) = sum((grad_shap(5:6)*inv_BK_t).*Normal,2);
  
  bK = RData.Vertices(1,:);
  BK = [RData.Vertices(2,:)-bK; ...
        RData.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  
  INR = zeros(1,3);
  dNR = zeros(1,3);
  switch(RData.EdgeLoc)
    case 1
      INR(1) = 0;
      INR(2) = dS/2;
      INR(3) = dS/2;
    case 2  
      INR(1) = dS/2;
      INR(2) = 0;
      INR(3) = dS/2;
    case 3
      INR(1) = dS/2;
      INR(2) = dS/2;
      INR(3) = 0;
  end
  dNR(1) = sum((grad_shap(1:2)*inv_BK_t).*Normal,2);
  dNR(2) = sum((grad_shap(3:4)*inv_BK_t).*Normal,2);
  dNR(3) = sum((grad_shap(5:6)*inv_BK_t).*Normal,2);
    
  % Compute entries of element penalty matrix
  
  if(LData.Element < RData.Element)
    gamma = 1;
  else
    gamma = -1;
  end
  
  Aloc(1,1) = gamma*( dNL(1)*INL(1)+s*INL(1)*dNL(1))/2;
  Aloc(1,3) = gamma*( dNL(1)*INL(2)+s*INL(1)*dNL(2))/2;
  Aloc(1,5) = gamma*( dNL(1)*INL(3)+s*INL(1)*dNL(3))/2;
  
  Aloc(1,7) = gamma*(-dNL(1)*INR(1)+s*INL(1)*dNR(1))/2;
  Aloc(1,9) = gamma*(-dNL(1)*INR(2)+s*INL(1)*dNR(2))/2;
  Aloc(1,11) = gamma*(-dNL(1)*INR(3)+s*INL(1)*dNR(3))/2;
  
  Aloc(2,2) = gamma*( dNL(1)*INL(1)+s*INL(1)*dNL(1))/2;
  Aloc(2,4) = gamma*( dNL(1)*INL(2)+s*INL(1)*dNL(2))/2;
  Aloc(2,6) = gamma*( dNL(1)*INL(3)+s*INL(1)*dNL(3))/2;
  
  Aloc(2,8) = gamma*(-dNL(1)*INR(1)+s*INL(1)*dNR(1))/2;
  Aloc(2,10) = gamma*(-dNL(1)*INR(2)+s*INL(1)*dNR(2))/2;
  Aloc(2,12) = gamma*(-dNL(1)*INR(3)+s*INL(1)*dNR(3))/2;
  
  Aloc(3,1) = gamma*( dNL(2)*INL(1)+s*INL(2)*dNL(1))/2;
  Aloc(3,3) = gamma*( dNL(2)*INL(2)+s*INL(2)*dNL(2))/2;
  Aloc(3,5) = gamma*( dNL(2)*INL(3)+s*INL(2)*dNL(3))/2;
  
  Aloc(3,7) = gamma*(-dNL(2)*INR(1)+s*INL(2)*dNR(1))/2;
  Aloc(3,9) = gamma*(-dNL(2)*INR(2)+s*INL(2)*dNR(2))/2;
  Aloc(3,11) = gamma*(-dNL(2)*INR(3)+s*INL(2)*dNR(3))/2;

  Aloc(4,2) = gamma*( dNL(2)*INL(1)+s*INL(2)*dNL(1))/2;
  Aloc(4,4) = gamma*( dNL(2)*INL(2)+s*INL(2)*dNL(2))/2;
  Aloc(4,6) = gamma*( dNL(2)*INL(3)+s*INL(2)*dNL(3))/2;
  
  Aloc(4,8) = gamma*(-dNL(2)*INR(1)+s*INL(2)*dNR(1))/2;
  Aloc(4,10) = gamma*(-dNL(2)*INR(2)+s*INL(2)*dNR(2))/2;
  Aloc(4,12) = gamma*(-dNL(2)*INR(3)+s*INL(2)*dNR(3))/2;
  
  Aloc(5,1) = gamma*( dNL(3)*INL(1)+s*INL(3)*dNL(1))/2;
  Aloc(5,3) = gamma*( dNL(3)*INL(2)+s*INL(3)*dNL(2))/2;
  Aloc(5,5) = gamma*( dNL(3)*INL(3)+s*INL(3)*dNL(3))/2;
  
  Aloc(5,7) = gamma*(-dNL(3)*INR(1)+s*INL(3)*dNR(1))/2;
  Aloc(5,9) = gamma*(-dNL(3)*INR(2)+s*INL(3)*dNR(2))/2;
  Aloc(5,11) = gamma*(-dNL(3)*INR(3)+s*INL(3)*dNR(3))/2;
  
  Aloc(6,2) = gamma*( dNL(3)*INL(1)+s*INL(3)*dNL(1))/2;
  Aloc(6,4) = gamma*( dNL(3)*INL(2)+s*INL(3)*dNL(2))/2;
  Aloc(6,6) = gamma*( dNL(3)*INL(3)+s*INL(3)*dNL(3))/2;
  
  Aloc(6,8) = gamma*(-dNL(3)*INR(1)+s*INL(3)*dNR(1))/2;
  Aloc(6,10) = gamma*(-dNL(3)*INR(2)+s*INL(3)*dNR(2))/2;
  Aloc(6,12) = gamma*(-dNL(3)*INR(3)+s*INL(3)*dNR(3))/2;
  
  Aloc(7,1) = gamma*( dNR(1)*INL(1)-s*INR(1)*dNL(1))/2;
  Aloc(7,3) = gamma*( dNR(1)*INL(2)-s*INR(1)*dNL(2))/2;
  Aloc(7,5) = gamma*( dNR(1)*INL(3)-s*INR(1)*dNL(3))/2;
  
  Aloc(7,7) = gamma*(-dNR(1)*INR(1)-s*INR(1)*dNR(1))/2;
  Aloc(7,9) = gamma*(-dNR(1)*INR(2)-s*INR(1)*dNR(2))/2;
  Aloc(7,11) = gamma*(-dNR(1)*INR(3)-s*INR(1)*dNR(3))/2;
  
  Aloc(8,2) = gamma*( dNR(1)*INL(1)-s*INR(1)*dNL(1))/2;
  Aloc(8,4) = gamma*( dNR(1)*INL(2)-s*INR(1)*dNL(2))/2;
  Aloc(8,6) = gamma*( dNR(1)*INL(3)-s*INR(1)*dNL(3))/2;
  
  Aloc(8,8) = gamma*(-dNR(1)*INR(1)-s*INR(1)*dNR(1))/2;
  Aloc(8,10) = gamma*(-dNR(1)*INR(2)-s*INR(1)*dNR(2))/2;
  Aloc(8,12) = gamma*(-dNR(1)*INR(3)-s*INR(1)*dNR(3))/2;
  
  Aloc(9,1) = gamma*( dNR(2)*INL(1)-s*INR(2)*dNL(1))/2;
  Aloc(9,3) = gamma*( dNR(2)*INL(2)-s*INR(2)*dNL(2))/2;
  Aloc(9,5) = gamma*( dNR(2)*INL(3)-s*INR(2)*dNL(3))/2;
  
  Aloc(9,7) = gamma*(-dNR(2)*INR(1)-s*INR(2)*dNR(1))/2; 
  Aloc(9,9) = gamma*(-dNR(2)*INR(2)-s*INR(2)*dNR(2))/2;
  Aloc(9,11) = gamma*(-dNR(2)*INR(3)-s*INR(2)*dNR(3))/2;
  
  Aloc(10,2) = gamma*( dNR(2)*INL(1)-s*INR(2)*dNL(1))/2;
  Aloc(10,4) = gamma*( dNR(2)*INL(2)-s*INR(2)*dNL(2))/2;
  Aloc(10,6) = gamma*( dNR(2)*INL(3)-s*INR(2)*dNL(3))/2;
  
  Aloc(10,8) = gamma*(-dNR(2)*INR(1)-s*INR(2)*dNR(1))/2; 
  Aloc(10,10) = gamma*(-dNR(2)*INR(2)-s*INR(2)*dNR(2))/2;
  Aloc(10,12) = gamma*(-dNR(2)*INR(3)-s*INR(2)*dNR(3))/2;
  
  Aloc(11,1) = gamma*( dNR(3)*INL(1)-s*INR(3)*dNL(1))/2;
  Aloc(11,3) = gamma*( dNR(3)*INL(2)-s*INR(3)*dNL(2))/2;
  Aloc(11,5) = gamma*( dNR(3)*INL(3)-s*INR(3)*dNL(3))/2;
  
  Aloc(11,7) = gamma*(-dNR(3)*INR(1)-s*INR(3)*dNR(1))/2;
  Aloc(11,9) = gamma*(-dNR(3)*INR(2)-s*INR(3)*dNR(2))/2;
  Aloc(11,11) = gamma*(-dNR(3)*INR(3)-s*INR(3)*dNR(3))/2;
  
  Aloc(12,2) = gamma*( dNR(3)*INL(1)-s*INR(3)*dNL(1))/2;
  Aloc(12,4) = gamma*( dNR(3)*INL(2)-s*INR(3)*dNL(2))/2;
  Aloc(12,6) = gamma*( dNR(3)*INL(3)-s*INR(3)*dNL(3))/2;
  
  Aloc(12,8) = gamma*(-dNR(3)*INR(1)-s*INR(3)*dNR(1))/2;
  Aloc(12,10) = gamma*(-dNR(3)*INR(2)-s*INR(3)*dNR(2))/2;
  Aloc(12,12) = gamma*(-dNR(3)*INR(3)-s*INR(3)*dNR(3))/2;
return