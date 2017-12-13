function Jloc = STIMA_InnPen_PDG(Edge,Normal,BdFlag,LData,RData,QuadRule, ...
                                 Shap,SHandle,varargin)
% STIMA_INNPEN_PDG Element stiffness matrix for interior penalty term.
%
%   JLOC = STIMA_INNPEN_PDG(EDGE,NORMAL,BDFLAG,LDATA,RDATA,QUARULE,SHAP, ...
%   SHANDLE) computes the entries of the element stiffness matrix for the
%   interior penalty term using the shape functions given by the function
%   handle SHAP.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   The integer BDFLAG denotes the boundary flag of the current edge. Note
%   that for interior edges only values larger than 0 are allowed.
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
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP is a function handle to the reference element shape functions.
%
%   SHANDLE is a function pointer to the edge weight function.
%
%   JLOC = STIMA_INNPEN_PDG(EDGE,NORMAL,BDFLAG,LDATA,RDATA,QUADRULE,SHAP, ...
%   SHANDLE,SPARAM) also handles the variable length argumet list SPARAM
%   to the function pointer SHANDLE.
%
%   Example:
%  
%   sigma = @(P1,P0,varargin)10/norm(P1-P0);
%   Jloc = STIMA_InnPen_PDG(Edge,Normal,BdFlags,LData,RData,gauleg(0,1,2), ...
%                           @shap_DGCR,sigma);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.x,1);  % Number of quadrature points
  
  % Compute jump weight
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  
  sigma = SHandle(P0,P1,varargin{:});
  
  % Compute edge length
  
  dS = norm(P1-P0);

  % Compute values of shape functions
     
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
  
  xr = zeros(nPts,2);
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
  
  % Evaluate shape functions
  
  NR = Shap(xr);
  NL = Shap(xl);
  nDofs = size(NL,2);
  
  % Preallocate memory
  
  Jloc = zeros(2*nDofs,2*nDofs);
  
  % Compute element penalty term
  
  for j1 = 1:nDofs
    for j2 = j1:nDofs
      Jloc(j1,j2) = sigma * sum(QuadRule.w.*NL(:,j1).*NL(:,j2))*dS;    
    end
  end
  for j1 = 1:nDofs
    for j2 = 1:nDofs
      Jloc(j1,j2+nDofs) = -sigma * sum(QuadRule.w.*NL(:,j1).*NR(:,j2))*dS;
    end
  end
  for j1 = 1:nDofs
    for j2 = j1:nDofs
      Jloc(nDofs+j1,nDofs+j2) = sigma * sum(QuadRule.w.*NR(:,j1).*NR(:,j2))*dS;    
    end
  end
      
  % Fill in lower trinagular part
    
  tri = triu(Jloc);
  Jloc = tril(transpose(tri),-1) + tri;
  
return
