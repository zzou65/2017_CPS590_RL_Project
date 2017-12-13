function L = assemNeu_QFE(Mesh,BdFlags,L,QuadRule,FHandle,varargin)
% ASSEMNEU_LFE Neumann boundary conditions.
%
%   L = ASSEMNEU_QFE(MESH,BDFLAG,L,QUADRULE,FHANDLE) incoporates the
%   Neumann boundary conditions with the data given by FHANDLE into the
%   right hand side load vector L. The boundary condition is only enforced
%   at edges whose BdFlag is equal to the integer BDFLAG. The 1D struct
%   QUADRULE is used to do the numerical integration along the edges.
%
%   L = ASSEMNEU_QFE(MESH,BDFLAG,L,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the boundary data function
%   FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first column
%                 specifies the left hand side element where the second column
%                 specifies the right hand side element.
%    EDGELOC      P-by-3 matrix connecting egdes to local edges of elements. 
%
%   Example:
%
%   L = assemNeu_QFE(Mesh,BdFlags,L,QuadRule,FHandle);
%
%   See also GET_BDEDGES, SHAP_QFE.

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nGauss = size(QuadRule.w,1);
  
  % Precompute shape functions
  
  N = shap_QFE([QuadRule.x zeros(nGauss,1)]);
  
  Lloc = zeros(3,1);
  for j1 = BdFlags
  
    % Extract Neumann edges
  
    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j1);
 
    for j2 = Loc'
    
      % Compute element map
      
      id_m = j2+nCoordinates;
      if(Mesh.Edge2Elem(j2,1))
        
        % Match orientation to left hand side element  
          
        Elem = Mesh.Edge2Elem(j2,1);
        EdgeLoc = Mesh.EdgeLoc(j2,1);
        id_s = Mesh.Elements(Elem,rem(EdgeLoc,3)+1);
        id_e = Mesh.Elements(Elem,rem(EdgeLoc+1,3)+1);
        
      else
        
        % Match orientation to right hand side element  
          
        Elem = Mesh.Edge2Elem(j2,2);
        EdgeLoc = Mesh.EdgeLoc(j2,2);
        id_s = Mesh.Elements(Elem,rem(EdgeLoc,3)+1);
        id_e = Mesh.Elements(Elem,rem(EdgeLoc+1,3)+1);
        
      end
      Q0 = Mesh.Coordinates(id_s,:);
      Q1 = Mesh.Coordinates(id_e,:);      
      x = ones(nGauss,1)*Q0+QuadRule.x*(Q1-Q0);
      dS = norm(Q1-Q0);
      
      % Evaluate Neumannn boundary data
      
      FVal = FHandle(x,j1,varargin{:});
      
      % Numerical integration along an edge
    
      Lloc(1) = sum(QuadRule.w.*FVal.*N(:,1))*dS;
      Lloc(2) = sum(QuadRule.w.*FVal.*N(:,2))*dS;
      Lloc(3) = sum(QuadRule.w.*FVal.*N(:,4))*dS;
      
      % Add contributions of Neumann data to load vector
      
      L(id_s) = L(id_s)+Lloc(1);
      L(id_m) = L(id_m)+Lloc(2);
      L(id_e) = L(id_e)+Lloc(3);
      
    end    
  end
  
return