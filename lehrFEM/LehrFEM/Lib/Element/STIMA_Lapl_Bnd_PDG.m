function Aloc = STIMA_Lapl_Bnd_PDG(Edge,Normal,BdFlag,Data,s,QuadRule, ...
                                   Shap,grad_Shap,varargin)
% STIMA_LAPL_BND_PDG Element stiffness matrix for boundary edge contributions.
%
%   ALOC = STIMA_LAPL_BND_PDG(EDGE,NORMAL,BDFLAG,DATA,S,QUADRULE,SHAP, ...
%   GRAD_SHAP) computes the entries of the element stiffness matrix for the
%   boundary edge contributions using the shape functions given by the
%   function handles SHAP and GRAD_SHAP.
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
%   The struct DATA contains the left or right hand side element data:
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
%   The integer S specifies wheter the diffusive fluxes are discretized in
%   a symmetric or anti-symmetric way:
%    +1 Antisymmetric discretization of diffusive fluxes
%    -1 Symmetric discretization of diffusive fluxes 
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP and GRAD_SHAP are function handles to the reference element shape
%   functions and their gradients.
%
%   Example:
%
%   Aloc = STIMA_InnPen_PDG(Edge,Normal,BdFlags,Data,1,gauleg(0,1,2), ...
%                           @shap_DGCR,@grad_shap_DGCR);
                    
%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
                          
  % Initialize constants

  nPts = size(QuadRule.x,1);  % Number of quadrature points

  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  % Compute edge length
  
  dS = norm(P1-P0);
  
  % Compute element mapping
  
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:)-bK; ...
        Data.Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  
  x = zeros(nPts,2);
  switch(Data.EdgeLoc)
    case 1
      if(Data.Match == +1)
        x(:,1) = 1 - QuadRule.x;
        x(:,2) = QuadRule.x;
      else
        x(:,1) = QuadRule.x;
        x(:,2) = 1 - QuadRule.x;
      end
    case 2
      if(Data.Match == +1)
        x(:,2) = 1 - QuadRule.x;  
      else
        x(:,2) = QuadRule.x;  
      end
    case 3
      if(Data.Match == +1)
        x(:,1) = QuadRule.x; 
      else
        x(:,1) = 1 - QuadRule.x;
      end;
  end
  
  % Evaluate shape functions and gradients
  
  N = Shap(x);
  grad_N = grad_Shap(x);
  nDofs = size(N,2);
  
  dN = zeros(nPts,nDofs);
  Normal = ones(nPts,1)*Normal;
  for j = 1:nDofs
    loc = 2*(j-1) + [1 2];
    dN(:,j) = sum((grad_N(:,loc)*inv_BK_t).*Normal,2);
  end
      
  % Preallocate memory
  
  Aloc = zeros(nDofs,nDofs);
  
  % Compute element stiffness matrix
    
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Aloc(j1,j2) = sum(QuadRule.w.*dN(:,j1).*N(:,j2))*dS + ...
                    s*sum(QuadRule.w.*N(:,j1).*dN(:,j2))*dS;
    end
  end
   
return
