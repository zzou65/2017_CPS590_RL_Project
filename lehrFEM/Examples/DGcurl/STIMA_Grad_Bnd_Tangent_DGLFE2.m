function Aloc = STIMA_Grad_Bnd_Tangent_DGLFE2(Edge,Normal,BdFlag,Data,s,varargin)
% STIMA_GRAD_BND_TANGENT_DGLFE2 Element stiffness matrix for boundary terms.
%
%   ALOC = STIMA_GRAD_BND_TANGENT_DGLFE2(EDGE,NORMAL,BDFLAG,DATA,S) computes the entries
%   of the element stiffness matrix for the boundary terms for vector valued function.
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
%   The structs DATA contains the left or right hand side element data:
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
  
  Aloc = zeros(6,6);
 
  % Compute normal derivatives and integrated shape functions
   
  dS = norm(Edge(2,:)-Edge(1,:));
 
  grad_shap = grad_shap_DGLFE([0 0]);
  
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:)-bK; ...
        Data.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  
  IN = zeros(1,6);
  dN = zeros(1,6);
  switch(Data.EdgeLoc)
    case 1
      IN(1) = 0;
      IN(2) = 0;
      IN(3) = dS/2;
      IN(4) = dS/2;
      IN(5) = dS/2;
      IN(6) = dS/2;
    case 2
      IN(1) = dS/2;
      IN(2) = dS/2;
      IN(3) = 0;
      IN(4) = 0;
      IN(5) = dS/2;
      IN(6) = dS/2;
    case 3
      IN(1) = dS/2;
      IN(2) = dS/2;
      IN(3) = dS/2;
      IN(4) = dS/2;
      IN(5) = 0;
      IN(6) = 0;
  end
  dN(1) = sum((grad_shap(1:2)*inv_BK_t).*Normal,2);
  dN(2) = sum((grad_shap(1:2)*inv_BK_t).*Normal,2);
  dN(3) = sum((grad_shap(3:4)*inv_BK_t).*Normal,2);
  dN(4) = sum((grad_shap(3:4)*inv_BK_t).*Normal,2);
  dN(5) = sum((grad_shap(5:6)*inv_BK_t).*Normal,2);
  dN(6) = sum((grad_shap(5:6)*inv_BK_t).*Normal,2);
  
  % Compute entries of element penalty matrix
  for i = 1:6
      for j = 1:6
          Aloc(i,j) = dN(i)*IN(j)+s*IN(i)*dN(j);
      end
  end
  
  for i = 1:2:6
      Aloc(i,:) = -Normal(2)*Aloc(i,:);
      Aloc(i+1,:) = Normal(1)*Aloc(i+1,:);
      Aloc(:,i) = -Normal(2)*Aloc(:,i);
      Aloc(:,i+1) = Normal(1)*Aloc(:,i+1); 
  end
return