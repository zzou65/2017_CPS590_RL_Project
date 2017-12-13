function Lloc = LOAD_Lapl_Bnd_PDG(Edge,Normal,BdFlag,Data,s,QuadRule,Shap, ...
                                  grad_Shap,SHandle,FHandle,varargin)
% LOAD_LAPL_BND_PDG Element load vector for boundary boundary load data.
%
%   LLOC = LOAD_LAPL_BND_PDG(EDGE,NORMAL,BDFLAG,DATA,S,QUADRULE,SHAP, ...
%   GRAD_SHAP,SHANDLE,FHANDLE) computes the entries of the element load
%   vector for the boundary load data using the shape functions given by
%   the function handle SHAP and GRAD_SHAP.
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
%   The integer S can specifies wheter the diffusive fluxes are discretized
%   in a symmetric or anti-symmetric way:
%    +1 Antisymmetric discretization of diffusive fluxes
%    -1 Symmetric discretization of diffusive fluxes
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP and GRAD_SHAP are function pointers to the reference element shape
%   functions and their gradients.
%
%   SHANDLE is a function pointer to the edge weight function.
%
%   FHANDLE is a function pointer to the load data.
%
%   LLOC = LOAD_LAPL_BND_PDG(EDGE,NORMAL,BDFLAG,DATA,S,QUADRULE,SHAP, ...
%   GRAD_SHAP,SHANDLE,FHANDLE,PARAM) also handles the variable length
%   argumet list PARAM to the function pointers SHANDLE and FHANDLE.
%
%   Example:
%
%   sigma = @(P0,P1,varargin)norm(P1-P0);
%   uD = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   Lloc = LOAD_Lapl_Bnd_PDG(Edge,Normal,BdFlags,Data,1,gauleg(0,1,2),@shap_DGCR, ...
%                            @grad_shap_DGCR,sigma,uD);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  
  nPts = size(QuadRule.x,1);
  
  % Compute value of jump weight function
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  sigma = SHandle(P0,P1,varargin{:});
 
  % Compute edge length
  
  dS = norm(P1-P0);
  
  % Compute element mapping
  
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:) - bK; ...
        Data.Vertices(3,:) - bK];
  inv_BK_t = transpose(inv(BK));
  
  xhat = zeros(nPts,2);
  switch(Data.EdgeLoc)
    case 1
      if(Data.Match == +1)
        xhat(:,1) = 1 - QuadRule.x;
        xhat(:,2) = QuadRule.x;
      else
        xhat(:,1) = QuadRule.x;
        xhat(:,2) = 1 - QuadRule.x;
      end
    case 2
      if(Data.Match == +1)
        xhat(:,2) = 1 - QuadRule.x;  
      else
        xhat(:,2) = QuadRule.x;  
      end
    case 3
      if(Data.Match == +1)
        xhat(:,1) = QuadRule.x; 
      else
        xhat(:,1) = 1 - QuadRule.x;
      end;
  end
  x = ones(nPts,1)*P0 + QuadRule.x*(P1-P0);
  
  % Evaluate shape functions and their gradients
  
  N = Shap(xhat); 
  grad_N = grad_Shap(xhat);
  
  nDofs = size(N,2);
  dN = zeros(nPts,nDofs);
  
  Normal = ones(nPts,1)*Normal;
  for i = 1:nDofs
    loc = 2*(i-1) + [1 2];
    dN(:,i) = sum((grad_N(:,loc)*inv_BK_t).*Normal,2);
  end
  
  % Compute function values
  
  FVal = FHandle(x,BdFlag,varargin{:});
  
  % Compute entries of element load vector
   
  Lloc = zeros(nDofs,1);
  for i = 1:nDofs
    Lloc(i) = sigma*sum(QuadRule.w.*FVal.*N(:,i))*dS + ...
              s*sum(QuadRule.w.*FVal.*dN(:,i))*dS;
  end
  
return