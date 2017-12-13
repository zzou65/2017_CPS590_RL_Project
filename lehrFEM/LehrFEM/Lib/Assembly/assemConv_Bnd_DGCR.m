function Bbnd = assemConv_Bnd_DGCR(Mesh,U,Lim,QuadRule,NumFlux,FHandle,varargin)
% ASSEMCONV_BND_DGCR Assemble non-linear boundary edge convection terms.
%
%   B = ASSEMCON_BND_DGCR(MESH,U,LIM,QUADRULE,NUMFLUX,FHANDLE) assembles
%   the global load vector from the local element contributions.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGES        N-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-3 matrix connecting egdes to local edges of elements. 
%    NORMALS      P-by-2 matrix specifying the normals on each edge. The
%                 normals on interior edges are chosen such that they point
%                 from the element with the lower number to the element
%                 with the higher number and on boundary edges such that
%                 they point outside the domain.
%    MATCH        P-by-2 matrix specifying wheter the edge orientation of
%                 the current edge matches the orientation of the left and
%                 right hand side element.
%
%   U denotes the FE ansatz function for which the boundary edge convection
%   term has to be assembled.
%
%   The integer LIM denotes the type of limiter to be used for the test
%   functions:
%    0  Piecwise constant test functions 
%    1  Piecewise linear test functions
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   The function handle NUMFLUX denotes the numerical flux function used
%   for the computation of the boundary edge convection term.
%    
%   The function handle FHANDLE denotes the Dirichlet boundary data of the
%   solution to the variational problem to be solved.
%
%   B = ASSEMCONV_BND_DGCR(MESH,U,LIM,QUADRULE,NUMFLUX,FHANDLE,PARAM) also
%   handles the variable length argument list PARAM to the numerical flux
%   function NUMFLUX and the boundary data FHANDLE.
%
%   Example:
%
%   uD = @(x,varargin)zeros(size(x,1),1);
%   B = assemConv_Inn_DGCR(Mesh,U,1,gauleg(0,1,2),@Upwind_LA,uD);
%
%   See also shap_DGCR.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements 
  nEdges = size(Mesh.Edges,1);        % Number of edges
  nPts = size(QuadRule.w,1);          % Number of quadrature points
  
  % Preallocate memory
  
  Bbnd = zeros(3*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Precompute shape function values
  
  u_shap = shap_DGCR([QuadRule.x zeros(nPts,1)]);
  switch(Lim)
    case 0
      v_shap = ones(size(u_shap))/3;  
    case 1
      v_shap = u_shap;  
  end
  
  % Assemble convection contributions
  
  N = zeros(nPts,3);
  for i = 1:nEdges
     
    if(BdFlags(i) < 0)
       
      % Compute element mapping  
        
      P0 = Mesh.Coordinates(Mesh.Edges(i,1),:);  
      P1 = Mesh.Coordinates(Mesh.Edges(i,2),:);
      Normal = Mesh.Normals(i,:);
              
      x = ones(nPts,1)*P0 + QuadRule.x*(P1-P0);
      
      dS = norm(P1-P0);
      
      % Extract left or right hand side element data
    
      if(Mesh.Edge2Elem(i,1) > 0)
        Elem = Mesh.Edge2Elem(i,1);
        Loc = Mesh.EdgeLoc(i,1);
        Match = Mesh.EdgeOrient(i,1);
      else
        Elem = Mesh.Edge2Elem(i,2);
        Loc = Mesh.EdgeLoc(i,2);
        Match = Mesh.EdgeOrient(i,2);
      end
      idx = 3*(Elem-1)+[1 2 3];
      
      % Compute value of DG solution and shape functions
      
      UVal = U(idx);
      switch(Loc)
        case 1
          if(Match == 1)
            LVal = UVal(1)*u_shap(:,3) + ...
                   UVal(2)*u_shap(:,1) + ...
                   UVal(3)*u_shap(:,2);  
            N(:,1) = v_shap(:,3);
            N(:,2) = v_shap(:,1);
            N(:,3) = v_shap(:,2);
          else
            LVal = UVal(1)*u_shap(:,3) + ...
                   UVal(2)*u_shap(:,2) + ...
                   UVal(3)*u_shap(:,1);            
            N(:,1) = v_shap(:,3);
            N(:,2) = v_shap(:,2);
            N(:,3) = v_shap(:,1);
          end
        case 2
          if(Match == 1)
            LVal = UVal(1)*u_shap(:,2) + ...
                   UVal(2)*u_shap(:,3) + ...
                   UVal(3)*u_shap(:,1);        
            N(:,1) = v_shap(:,2);
            N(:,2) = v_shap(:,3);
            N(:,3) = v_shap(:,1);
          else
            LVal = UVal(1)*u_shap(:,1) + ...
                   UVal(2)*u_shap(:,3) + ...
                   UVal(3)*u_shap(:,2); 
            N(:,1) = v_shap(:,1);
            N(:,2) = v_shap(:,3);
            N(:,3) = v_shap(:,2);
          end
        case 3
          if(Match == 1)
            LVal = UVal(1)*u_shap(:,1) + ...
                   UVal(2)*u_shap(:,2) + ...
                   UVal(3)*u_shap(:,3);
            N(:,1) = v_shap(:,1);
            N(:,2) = v_shap(:,2);
            N(:,3) = v_shap(:,3);
          else
            LVal = UVal(1)*u_shap(:,2) + ...
                   UVal(2)*u_shap(:,1) + ...
                   UVal(3)*u_shap(:,3);        
            N(:,1) = v_shap(:,2);
            N(:,2) = v_shap(:,1);
            N(:,3) = v_shap(:,3);
          end
      end      
      
      % Compute value of boundary data
    
      RVal = FHandle(x,BdFlags(i),varargin{:});
    
      % Evaluate numerical flux function
      
      HVal = NumFlux(x,LVal,RVal,Normal,varargin{:});
      
      % Add element contributions to stiffness matrix
    
      Bbnd(idx(1)) = Bbnd(idx(1)) + sum(QuadRule.w.*HVal.*N(:,1))*dS;
      Bbnd(idx(2)) = Bbnd(idx(2)) + sum(QuadRule.w.*HVal.*N(:,2))*dS;
      Bbnd(idx(3)) = Bbnd(idx(3)) + sum(QuadRule.w.*HVal.*N(:,3))*dS;
   
    end  
  end
  
return
