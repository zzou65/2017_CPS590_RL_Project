function L = assemLoad_Bnd_PDG2_vec(Mesh,BdFlag,dim,EHandle,varargin)
%ASSEMLOAD_BND_PDG2_VEC Assemble vector-valued DG boundary edge contributions
%
%   L = ASSEMLOAD_BND_PDG2_VEC(MESH,BDFLAG,DIM,EHANDLE) assembles the
%   global load vector from local edge contributions given by the function 
%   handle EHANDLE.  This code is for vector-valued functions, for example
%   mixed DG, with the same basis functions in each component of the image.
%   The structure of the discrete vectors should be such that all
%   coefficients corresponding to the first dimension come first, then the 
%   second dimension, and so on.
%
%   L = ASSEMLOAD_BND_PDG2_VEC(MESH,BDFLAG,DIM,EHANDLE,EPARAM) passes the
%   variable-length argument list EPARAM to the function handle EHANDLE
%   during the assembly process.
%
%   BDFLAG is a scalar that specifies which part of the boundary is to be
%   processed.  The edges with flags equal to BDFLAG are taken into
%   account.  If BDFLAG is empty, all edges with negative boundary flags
%   are processed.
%
%   DIM is the dimension of the image space.  For example, for mixed DG,
%   the image space consists of one dimension for the function value and
%   two for the gradient, amounting to a total of three.
%
%   EHANDLE must take at least the arguments:
%     EDGE        Coordinates of edge
%     NORMAL      Normal vector
%     EDGEDATA    Struct containing arbitrary edge information
%     DATA        Struct containing information on adjacent element
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the 
%                 mesh.
%    EDGES        L-by-2 matrix specifying all edges of the mesh.
%    BDFLAGS      L-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies whether the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    L-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      L-by-3 or P-by-4 matrix connecting egdes to local edges
%                 of elements. 
%    NORMALS      L-by-2 matrix specifying the normals on each edge. The
%                 normals on interior edges are chosen such that they point
%                 from the element with the lower number to the element
%                 with the higher number and on boundary edges such that
%                 they point outside the domain.
%    MATCH        L-by-2 matrix specifying whether the edge orientation of
%                 the current edge matches the orientation of the left and
%                 right hand side element.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%    EDGEDATA     L-by-1 structure array containing at least the fields:
%       NDOFS       The total number of degrees of freedom on the two
%                   elements adjacent to a given edge.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, 

  % Initialize constants
  nEdges = size(Mesh.Edges,1);      % Number of edges
  
  % Find relevant edges
  if(isempty(BdFlag))
    doEdge = Mesh.BdFlags<0;
  else
    doEdge = Mesh.BdFlags==BdFlag;
  end
  
  % Count degrees of freedom on all elements and calculate the number of
  % matrix entries.
  nDofsSum = cumsum([0,[Mesh.ElemData.nDofs]]); % Number of degrees of freedom on all preceding elements
  nDofsTot = nDofsSum(end);             % Number of degrees of freedom per dimension
  numel = dim*nDofsTot;                 % Number of vector entries
  
  % Preallocate memory
  L = zeros(numel,1);
  
  % Assemble element contributions
  dimOffset = nDofsTot*(0:dim-1);
  for i = 1:nEdges
    
    % Check for boundary edge  
      
    if(doEdge(i))
      
      % Extract edge data
      Edge = Mesh.Coordinates(Mesh.Edges(i,:),:);
      Normal = Mesh.Normals(i,:);

      % Extract left or right hand side element data
      if(Mesh.Edge2Elem(i,1) > 0)
        Data.Element = Mesh.Edge2Elem(i,1);
        Data.ElemData = Mesh.ElemData(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,1);
        Data.Match = Mesh.EdgeOrient(i,1);
      else
        Data.Element = Mesh.Edge2Elem(i,2);
        Data.ElemData = Mesh.ElemData(Data.Element);
        Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
        Data.EdgeLoc = Mesh.EdgeLoc(i,2);
        Data.Match = Mesh.EdgeOrient(i,2);
      end
            
      % Compute element contributions
      Lloc = EHandle(Edge,Normal,Mesh.EdgeData(i), ...
                     Data,varargin{:}); 
    
      % Add element contributions to load vector
      idx0 = nDofsSum(Data.Element) + (1:Data.ElemData.nDofs);
      idx = zeros(1,dim*Data.ElemData.nDofs);
      idx(:) = idx0(ones(dim,1),:)' + dimOffset(ones(Data.ElemData.nDofs,1),:);
      L(idx) = L(idx) + Lloc;      
          
    end
  end

return