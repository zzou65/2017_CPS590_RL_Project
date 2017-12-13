function L = assemLoad_Bnd_DG(Mesh,EHandle,varargin)
% ASSEMLOAD_BND_DG Assemble discontinuous Crouzeix-Raviart FE boundary edge
% contributions.
%
%   A = ASSEMLOAD_BND_DG(MESH,EHANDLE) assembles the global load vector
%   from the local element contributions given by the function handle
%   EHANDLE.
%
%   A = ASSEMLOAD_BND_DG(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
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
%   Example:
%
%   sigma = @(P0,P1,varargin)10/norm(P1-P0);
%   uD = @(x,varargin)zeros(size(x,1),1);
%   Lbnd = assemLoad_Bnd_DG(Mesh,@LOAD_Bnd_DGCR,gauleg(0,1,2),0,sigma,uD);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nEdges = size(Mesh.Edges,1);        % Number of edges
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Preallocate memory
  
  L = zeros(3*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Assemble element contributions
  
  for i = 1:nEdges
    
    % Check for boundary edge  
      
    if(BdFlags(i) < 0)
        
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Normal = Mesh.Normals(i,:);

      % Extract left or right hand side element data
    
      if(Mesh.Edge2Elem(i,1) > 0)
        Data.Element = Mesh.Edge2Elem(i,1);
        Data.ElemFlag = ElemFlag(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,1);
        Data.Match = Mesh.EdgeOrient(i,1);
      else
        Data.Element = Mesh.Edge2Elem(i,2);
        Data.ElemFlag = ElemFlag(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,2);
        Data.Match = Mesh.EdgeOrient(i,2);
      end
            
      % Compute element contributions
    
      Lloc = EHandle(Edge,Normal,BdFlags(i), ...
                     Data,varargin{:});
    
      % Add element contributions to load vector
    
      idx = 3*(Data.Element-1) + [1 2 3];
    
      L(idx) = L(idx) + Lloc;      
          
    end
  end

return