function Aloc = STIMA_Bnd_DGLFE(Edge,Normal,BdFlag,Data,s,varargin)
% STIMA_BND_DGLFE Element stiffness matrix for boundary terms.
%
%   ALOC = STIMA_BND_DGLFE(EDGE,NORMAL,BDFLAG,DATA,S) computes the entries
%   of the element stiffness matrix for the boundary terms.
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
%
%   Example:
%
%   Aloc = STIMA_Bnd_DGLFE(Edge,Normal,BdFlags,Data,1);
%
%   See also grad_shap_DGLFE.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(3,3);
 
  % Compute normal derivatives and integrated shape functions
   
  dS = norm(Edge(2,:)-Edge(1,:));
 
  grad_shap = grad_shap_DGLFE([0 0]);
  
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:)-bK; ...
        Data.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  
  IN = zeros(1,3);
  dN = zeros(1,3);
  switch(Data.EdgeLoc)
    case 1
      IN(1) = 0;
      IN(2) = dS/2;
      IN(3) = dS/2;
    case 2
      IN(1) = dS/2;
      IN(2) = 0;
      IN(3) = dS/2;
    case 3
      IN(:,1) = dS/2;
      IN(:,2) = dS/2;
      IN(:,3) = 0;
  end
  dN(1) = sum((grad_shap(1:2)*inv_BK_t).*Normal,2);
  dN(2) = sum((grad_shap(3:4)*inv_BK_t).*Normal,2);
  dN(3) = sum((grad_shap(5:6)*inv_BK_t).*Normal,2);
  
  % Compute entries of element penalty matrix
  
  Aloc(1,1) = dN(1)*IN(1)+s*IN(1)*dN(1);
  Aloc(1,2) = dN(1)*IN(2)+s*IN(1)*dN(2);
  Aloc(1,3) = dN(1)*IN(3)+s*IN(1)*dN(3);
  
  Aloc(2,1) = dN(2)*IN(1)+s*IN(2)*dN(1);
  Aloc(2,2) = dN(2)*IN(2)+s*IN(2)*dN(2);
  Aloc(2,3) = dN(2)*IN(3)+s*IN(2)*dN(3);
    
  Aloc(3,1) = dN(3)*IN(1)+s*IN(3)*dN(1);
  Aloc(3,2) = dN(3)*IN(2)+s*IN(3)*dN(2);
  Aloc(3,3) = dN(3)*IN(3)+s*IN(3)*dN(3);
        
return
