function varargout = assemMat_Bnd_LieW1F(Mesh,BdFlagsParam,EHandle,varargin)
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
%   See also set_Rows, set_Cols.

%   Copyright 2006-2009  Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

nElements = size(Mesh.Elements,1);  % Number of elements in the mesh
nEdges = size(Mesh.Edges,1);        % Number of edges in the mesh

% Allocate memory

I = zeros(9*nEdges,1);
J = zeros(9*nEdges,1);
A = zeros(9*nEdges,1);

% Check for element flags

if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag;
else
    ElemFlag = zeros(nElements,1);
end
BdFlags = Mesh.BdFlags;

% Assemble element contributions

loc = 1:9;
last = 0;

for j1=BdFlagsParam
    
    % Extract certain edges
    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j1);

    for i = Loc'

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

        Aloc = EHandle(Edge,Normal,BdFlags(i), ...
            Data,varargin{:});

        % Add element contributions to stiffness matrix
        idx = [Mesh.Vert2Edge(Mesh.Elements(Data.Element,2),Mesh.Elements(Data.Element,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(Data.Element,3),Mesh.Elements(Data.Element,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(Data.Element,1),Mesh.Elements(Data.Element,2))];

        % Determine the orientation
        if(Mesh.Edges(idx(1),1)==Mesh.Elements(Data.Element,2)),  p1 = 1;  else    p1 = -1;  end
        if(Mesh.Edges(idx(2),1)==Mesh.Elements(Data.Element,3)),  p2 = 1;  else    p2 = -1;  end
        if(Mesh.Edges(idx(3),1)==Mesh.Elements(Data.Element,1)),  p3 = 1;  else    p3 = -1;  end

        Peori = diag([p1 p2 p3]); % scaling matrix taking into account orientations
        Aloc = Peori*Aloc*Peori;

        I(loc) = set_Rows(idx,3);
        J(loc) = set_Cols(idx,3);
        A(loc) = Aloc(:);

        loc = loc + 9;
        last = last + 9;

    end

end

% Assign output arguments

loc = 1:last;
if(nargout > 1)
    varargout{1} = I(loc);
    varargout{2} = J(loc);
    varargout{3} = A(loc);
else
    varargout{1} = sparse(I(loc),J(loc),A(loc),nEdges,nEdges);
end

return
