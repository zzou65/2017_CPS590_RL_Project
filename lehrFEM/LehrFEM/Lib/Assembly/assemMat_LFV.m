function varargout = assemMat_LFV(Mesh,EHandle,varargin)
% ASSEMMAT_LFV Assemble linear FV contributions.
%
%   A = ASSEMMAT_LFV(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_LFV(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process.
%
%   [I,J,A] = ASSEMMAT_LFV(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGES        L-by-2 matrix specifying the edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix specifying the edge connecting vertices.
%    BDFLAGS	  L-by-1 matrix specifying boundary flags for the edges.
%    MIDPOINTS    L-by-2 matrix specifying the border points on the edges.
%    CENTERPOINTS N-by-2 matrix specifying the center points of the elements.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   Mesh = add_Edges(Mesh);
%   Mesh = add_MidPoints(Mesh);
%   EHandle = @STIMA_GenLapl_LFV;
%   A = assemMat_LFV(Mesh,EHandle,@(x)1);
%  
%   See also set_Rows, set_Cols.

%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  I = zeros(9*nElements,1);
  J = zeros(9*nElements,1);
  A = zeros(9*nElements,1);
  
  % Check for element flags
  
  %if(isfield(Mesh,'ElemFlag')),
  %  flags = Mesh.ElemFlag; 
  %else
  %  flags = zeros(nElements,1);
  %end

  % Check for boundary flags
  if isfield(Mesh,'BdFlags')
    bdFlags = Mesh.BdFlags;
  else
    bdFlags = zeros(size(Mesh.Edges,1),1);
  end

  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices, midpoints, centerpoint
    % and boundary flags of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
    MidPoints(1,:) = Mesh.MidPoints(Mesh.Vert2Edge(idx(1),idx(2)),:);
    MidPoints(2,:) = Mesh.MidPoints(Mesh.Vert2Edge(idx(2),idx(3)),:);
    MidPoints(3,:) = Mesh.MidPoints(Mesh.Vert2Edge(idx(3),idx(1)),:);
    Center = Mesh.CenterPoints(i,:);

    bF(1) = bdFlags(Mesh.Vert2Edge(idx(1),idx(2)));
    bF(2) = bdFlags(Mesh.Vert2Edge(idx(2),idx(3)));
    bF(3) = bdFlags(Mesh.Vert2Edge(idx(3),idx(1)));
   
    % Compute element contributions
    
    Aloc = EHandle(Vertices,MidPoints,Center,Mesh.Type,bF,varargin{:});
    
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(idx,3);
    J(loc) = set_Cols(idx,3);
    A(loc) = Aloc(:);
    loc = loc+9;
    
  end
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);      
  end
  
return
