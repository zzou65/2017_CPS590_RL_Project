function Binn = assemConv_Inn_DGCR(Mesh,U,Lim,QuadRule,NumFlux,varargin)
% ASSEMCONV_INN_DGCR Assemble non-linear interior edge convection terms.
%
%   B = ASSEMCON_INN_DGCR(MESH,U,LIM,QUADRULE,NUMFLUX) assembles the global
%   load vector from the local element contributions.
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
%   U deotes the FE ansatz function for which the interior edge convection
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
%   for the computation of the interior edge convection term.
%    
%   B = ASSEMCONV_INN_DGCR(MESH,U,LIM,QUADRULE,NUMFLUX,PARAM) also handles
%   the variable length argument list PARAM to the numerical flux function
%   NUMFLUX.
%
%   Example:
%
%   B = assemConv_Inn_DGCR(Mesh,U,1,gauleg(0,1,2),@Upwind_LA);
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
  
  Binn = zeros(3*nElements,1);
  
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
 
  NL = zeros(nPts,3);
  NR = zeros(nPts,3);
  for i = 1:nEdges
     
    if(BdFlags(i) >= 0)
       
      % Compute element mapping  
        
      P0 = Mesh.Coordinates(Mesh.Edges(i,1),:);  
      P1 = Mesh.Coordinates(Mesh.Edges(i,2),:);
      Normal = Mesh.Normals(i,:);
              
      x = ones(nPts,1)*P0 + QuadRule.x*(P1-P0);
      
      dS = norm(P1-P0);
            
      % Compute value of DG solution on left and right hand side element
      
      if(Mesh.Edge2Elem(i,1) < Mesh.Edge2Elem(i,2))
        LElem = Mesh.Edge2Elem(i,1);
        LLoc = Mesh.EdgeLoc(i,1);
        LMatch = Mesh.EdgeOrient(i,1);
        RElem = Mesh.Edge2Elem(i,2);
        RLoc = Mesh.EdgeLoc(i,2);
        RMatch = Mesh.EdgeOrient(i,2);
      else
        LElem = Mesh.Edge2Elem(i,2);
        LLoc = Mesh.EdgeLoc(i,2);
        LMatch = Mesh.EdgeOrient(i,2);
        RElem = Mesh.Edge2Elem(i,1);
        RLoc = Mesh.EdgeLoc(i,1);
        RMatch = Mesh.EdgeOrient(i,1);
      end
      Lidx = 3*(LElem-1)+[1 2 3];
      UL = U(Lidx);
      Ridx = 3*(RElem-1)+[1 2 3];
      UR = U(Ridx);

      % Compute traces of DG solution and shape functions
      
      switch(LLoc)
        case 1
          if(LMatch == 1)     
            LVal = UL(1)*u_shap(:,3) + ...
                   UL(2)*u_shap(:,1) + ...
                   UL(3)*u_shap(:,2);  
            NL(:,1) = v_shap(:,3);
            NL(:,2) = v_shap(:,1);
            NL(:,3) = v_shap(:,2);
          else
            LVal = UL(1)*u_shap(:,3) + ...
                   UL(2)*u_shap(:,2) + ...
                   UL(3)*u_shap(:,1);               
            NL(:,1) = v_shap(:,3);
            NL(:,2) = v_shap(:,2);
            NL(:,3) = v_shap(:,1);
          end
        case 2
          if(LMatch == 1)
            LVal = UL(1)*u_shap(:,2) + ...
                   UL(2)*u_shap(:,3) + ...
                   UL(3)*u_shap(:,1);        
            NL(:,1) = v_shap(:,2);
            NL(:,2) = v_shap(:,3);
            NL(:,3) = v_shap(:,1);
          else
            LVal = UL(1)*u_shap(:,1) + ...
                   UL(2)*u_shap(:,3) + ...
                   UL(3)*u_shap(:,2);       
            NL(:,1) = v_shap(:,1);
            NL(:,2) = v_shap(:,3);
            NL(:,3) = v_shap(:,2);
          end
        case 3
          if(LMatch == 1)
            LVal = UL(1)*u_shap(:,1) + ...
                   UL(2)*u_shap(:,2) + ...
                   UL(3)*u_shap(:,3);  
            NL(:,1) = v_shap(:,1);
            NL(:,2) = v_shap(:,2);
            NL(:,3) = v_shap(:,3);
          else
            LVal = UL(1)*u_shap(:,2) + ...
                   UL(2)*u_shap(:,1) + ...
                   UL(3)*u_shap(:,3);        
            NL(:,1) = v_shap(:,2);
            NL(:,2) = v_shap(:,1);
            NL(:,3) = v_shap(:,3);
          end
      end      
      
      switch(RLoc)
        case 1
          if(RMatch == 1)
            RVal = UR(1)*u_shap(:,3) + ...
                   UR(2)*u_shap(:,1) + ...
                   UR(3)*u_shap(:,2);
            NR(:,1) = v_shap(:,3);
            NR(:,2) = v_shap(:,1);
            NR(:,3) = v_shap(:,2);
          else
            RVal = UR(1)*u_shap(:,3) + ...
                   UR(2)*u_shap(:,2) + ...
                   UR(3)*u_shap(:,1);       
            NR(:,1) = v_shap(:,3);
            NR(:,2) = v_shap(:,2);
            NR(:,3) = v_shap(:,1);
          end
        case 2
          if(RMatch == 1)
            RVal = UR(1)*u_shap(:,2) + ...
                   UR(2)*u_shap(:,3) + ...
                   UR(3)*u_shap(:,1);         
            NR(:,1) = v_shap(:,2);
            NR(:,2) = v_shap(:,3);
            NR(:,3) = v_shap(:,1);
          else
            RVal = UR(1)*u_shap(:,1) + ...
                   UR(2)*u_shap(:,3) + ...
                   UR(3)*u_shap(:,2);       
            NR(:,1) = v_shap(:,1);
            NR(:,2) = v_shap(:,3);
            NR(:,3) = v_shap(:,2);
          end
        case 3
          if(RMatch == 1)
            RVal = UR(1)*u_shap(:,1) + ...
                   UR(2)*u_shap(:,2) + ...
                   UR(3)*u_shap(:,3);  
            NR(:,1) = v_shap(:,1);
            NR(:,2) = v_shap(:,2);
            NR(:,3) = v_shap(:,3);
          else
            RVal = UR(1)*u_shap(:,2) + ...
                   UR(2)*u_shap(:,1) + ...
                   UR(3)*u_shap(:,3);        
            NR(:,1) = v_shap(:,2);
            NR(:,2) = v_shap(:,1);
            NR(:,3) = v_shap(:,3);
          end
      end      
            
      % Evaluate numerical flux function
      
      HVal = NumFlux(x,LVal,RVal,Normal,varargin{:});
      
      % Add convection contributions
      
      Binn(Lidx(1)) = Binn(Lidx(1)) + sum(QuadRule.w.*HVal.*NL(:,1))*dS;
      Binn(Lidx(2)) = Binn(Lidx(2)) + sum(QuadRule.w.*HVal.*NL(:,2))*dS;
      Binn(Lidx(3)) = Binn(Lidx(3)) + sum(QuadRule.w.*HVal.*NL(:,3))*dS;
      
      Binn(Ridx(1)) = Binn(Ridx(1)) - sum(QuadRule.w.*HVal.*NR(:,1))*dS;
      Binn(Ridx(2)) = Binn(Ridx(2)) - sum(QuadRule.w.*HVal.*NR(:,2))*dS;
      Binn(Ridx(3)) = Binn(Ridx(3)) - sum(QuadRule.w.*HVal.*NR(:,3))*dS;
         
    end
      
  end

return
