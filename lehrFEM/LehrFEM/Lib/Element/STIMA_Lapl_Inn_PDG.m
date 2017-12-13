function Aloc = STIMA_Lapl_Inn_PDG(Edge,Normal,BdFlag,LData,RData,s, ...
                                   QuadRule,Shap,grad_Shap,varargin)
% STIMA_LAPL_INN_PDG Element stiffness matrix for interior edge contributions.
%
%   ALOC = STIMA_LAPL_INN_PDG(EDGE,NORMAL,BDFLAG,LDATA,RDATA,S,QUARULE, ... 
%   SHAP,GRAD_SHAP) computes the entries of the element stiffness matrix for
%   the interior edge contributions using the shape functions given by the
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
%   Aloc = STIMA_Lapl_Inn_PDG(Edge,Normal,BdFlags,LData,RData,1, ...
%                             gauleg(0,1,2),@shap_DGCR,@grad_shap_DGCR);

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
     
  bK_L = LData.Vertices(1,:);
  BK_L = [LData.Vertices(2,:)-bK_L; ...
          LData.Vertices(3,:)-bK_L];
  inv_BK_t_L = transpose(inv(BK_L));
  
  xl = zeros(nPts,2);
  switch(LData.EdgeLoc)
    case 1
      if(LData.Match == +1)
        xl(:,1) = 1 - QuadRule.x;
        xl(:,2) = QuadRule.x;
      else
        xl(:,1) = QuadRule.x;
        xl(:,2) = 1 - QuadRule.x;
      end
    case 2
      if(LData.Match == +1)
        xl(:,2) = 1 - QuadRule.x;  
      else
        xl(:,2) = QuadRule.x;  
      end
    case 3
      if(LData.Match == +1)
        xl(:,1) = QuadRule.x;  
      else
        xl(:,1) = 1 - QuadRule.x;  
      end 
  end
  
  bK_R = RData.Vertices(1,:);
  BK_R = [RData.Vertices(2,:)-bK_R; ...
          RData.Vertices(3,:)-bK_R];
  inv_BK_t_R = transpose(inv(BK_R));
  
  xr = zeros(nPts);
  switch(RData.EdgeLoc)     
    case 1
      if(RData.Match == +1)
        xr(:,1) = 1 - QuadRule.x;
        xr(:,2) = QuadRule.x; 
      else
        xr(:,1) = QuadRule.x;
        xr(:,2) = 1 - QuadRule.x;  
      end
    case 2
       if(RData.Match == +1)
         xr(:,2) = 1 - QuadRule.x;  
       else
         xr(:,2) = QuadRule.x;  
       end
    case 3
       if(RData.Match == +1)
         xr(:,1) = QuadRule.x;  
       else
         xr(:,1) = 1 - QuadRule.x;  
       end
  end
 
  % Evaluate shape functions and gradients
  
  NL = Shap(xl);
  grad_NL = grad_Shap(xl);
  NR = Shap(xr);
  grad_NR = grad_Shap(xr);
  nDofs = size(NR,2);
  
  dNL = zeros(nPts,nDofs);
  dNR = zeros(nPts,nDofs);
  Normal = ones(nPts,1)*Normal;
  for j = 1:nDofs
    loc = 2*(j-1) + [1 2];
    dNL(:,j) = sum((grad_NL(:,loc)*inv_BK_t_L).*Normal,2);
    dNR(:,j) = sum((grad_NR(:,loc)*inv_BK_t_R).*Normal,2);
  end

  % Preallocate memory
  
  Aloc = zeros(2*nDofs,2*nDofs);
  
  % Compute element penalty term
  
  if(RData.Element > LData.Element)
    gamma = +1/2;
  else
    gamma = -1/2;  
  end
     
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Aloc(j1,j2) = + gamma*sum(QuadRule.w.*dNL(:,j1).*NL(:,j2))*dS ...
                    + s*gamma*sum(QuadRule.w.*NL(:,j1).*dNL(:,j2))*dS;
    end
  end
      
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Aloc(j1+nDofs,j2) = + gamma*sum(QuadRule.w.*dNR(:,j1).*NL(:,j2))*dS ...
                          - s*gamma*sum(QuadRule.w.*NR(:,j1).*dNL(:,j2))*dS;      
    end
  end
  
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Aloc(j1,j2+nDofs) = - gamma*sum(QuadRule.w.*dNL(:,j1).*NR(:,j2))*dS ...
                          + s*gamma*sum(QuadRule.w.*NL(:,j1).*dNR(:,j2))*dS;
    end                               
  end
  
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Aloc(j1+nDofs,j2+nDofs) = - gamma*sum(QuadRule.w.*dNR(:,j1).*NR(:,j2))*dS ...
                                - s*gamma*sum(QuadRule.w.*NR(:,j1).*dNR(:,j2))*dS;  
    end
  end
  
return