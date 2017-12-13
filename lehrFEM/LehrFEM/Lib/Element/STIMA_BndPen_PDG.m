function Jloc = STIMA_BndPen_PDG(Edge,Normal,BdFlag,Data,QuadRule,Shap, ...
                                 SHandle,varargin)
% STIMA_BNDPEN_PDG Element stiffness matrix for boundary penalty term.
%
%   JLOC = STIMA_BNDPEN_PDG(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,SHAP,SHANDLE)
%   computes the entries of the element stiffness matrix for the boundary
%   penalty term using the shape functions given by the function handle
%   SHAP.
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
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP is a function handle to the reference element shape functions.
%
%   SHANDLE is a function pointer to the edge weight function.
%
%   JLOC = STIMA_BNDPEN_PDG(EDGE,NORMAL,BDFLAG,LDATA,RDATA,QUADRULE,SHAP, ...
%   SHANDLE,SPARAM) also handles the variable length argumet list SPARAM to
%   the function pointer SHANDLE.
%
%   Example:
%
%   sigma = @(P1,P0,varargin)10/norm(P1-P0);
%   Jloc = STIMA_InnPen_PDG(Edge,Normal,BdFlags,Data,gauleg(0,1,2), ...
%                           @shap_DGCR,sigma);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nPts = size(QuadRule.x,1);  % Number of qudrature points

  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  sigma = SHandle(P0,P1,varargin{:});
  
  % Compute edge length
  
  dS = norm(P1-P0);
  
  % Compute element mapping
 
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
  
  % Evaluate shape functions
  
  N = Shap(x);
  nDofs = size(N,2);
  
  % Preallocate memory
  
  Jloc = zeros(nDofs,nDofs);
  
  % Compute element penalty term
  
  for j1 = 1:nDofs
    for j2 = j1:nDofs
      Jloc(j1,j2) = sigma*sum(QuadRule.w.*N(:,j1).*N(:,j2))*dS;
    end
  end
  
  % Fill in lower trinagular part
    
  tri = triu(Jloc);
  Jloc = tril(transpose(tri),-1) + tri;
  
return
