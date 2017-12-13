function L = assemLoad_PML_Bnd(Mesh,BdFlags,L,QuadRule,FHandle, eps_Handle,varargin)
% ASSEMLOAD_PML_Bnd PML righthandside .
%
%   L = ASSEMLOAD_PML(MESH,BDFLAG,L,QUADRULE,FHANDLE, eps_Handle) incoporates the
%   Neumann-like boundary term with the data given by FHANDLE into the
%   right hand side load vector L. The boundary condition is only enforced
%   at edges whose BdFlag is equal to the integer BDFLAG. The 1D struct
%   QUADRULE is used to do the numerical integration along the edges.
%
%   L = ASSEMLoad_PML_Bnd(MESH,BDFLAG,L,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the boundary data function
%   FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first column
%                 specifies the left hand side element where the second column
%                 specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of elements. 
%
%   Example:
%
%   L = assemLoad_PML_Bnd(Mesh,BdFlags,L,QuadRule,FHandle, epsHandle);
%
%   See also GET_BDEDGES, SHAP_LFE.

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia & Holger Heumann &
%   Nana
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nGauss = size(QuadRule.w,1);
  
  % Precompute shape functions
  
  N = shap_LFE([QuadRule.x zeros(nGauss,1)]);
  
  Lloc = zeros(2,1);
  for j1 = BdFlags
  
    % Extract Neumann edges
  
    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j1);
    
    for j2 = Loc'
        
      % Compute element map
      
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
      tang=(Q1-Q0);
      dS = norm(tang);
      Normal=ones(size(x,2),1)*[-tang(2) tang(1)];
      
      % Evaluate boundary data, evaluate Matrix products
      
      FVal = FHandle(x);
      epsVal=EpsHandle(x);
      Val=epsVal(:,1).*Fval(:,1).*Normal(:,1)+...
          epsVal(:,2).*Fval(:,2).*Normal(:,1)+...
          epsVal(:,3).*Fval(:,1).*Normal(:,2)+...
          epsVal(:,4).*Fval(:,2).*Normal(:,2);
      
      % Numerical integration along an edge
    
      Lloc(1) = sum(QuadRule.w.*Val.*N(:,1));
      Lloc(2) = sum(QuadRule.w.*Val.*N(:,2));
  
      % Add contributions of Neumann data to load vector
      
      L(id_s) = L(id_s)+Lloc(1);
      L(id_e) = L(id_e)+Lloc(2);
      
    end    
  end
  
return