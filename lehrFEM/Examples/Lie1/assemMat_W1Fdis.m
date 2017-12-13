function varargout = assemMat_W1Fdis(Mesh,Sigma,varargin)
% ASSEMMAT_W1F Assembly for *edge elements* in 2D using sigma as piecewise 
%  constant ( on cells) parameter
%
%   A = ASSEMMAT_W1F(MESH,Sigma) assembles the global matrix from the
%   local element contributions given by vertice array Sigma and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMAT_W1F(MESH,EHANDLE) assembles the global matrix from
%   the local element contributions given by the function handle EHANDLE
%   and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    VERT2EDGE    Edge numbers associated with pairs of vertices 
%                 (sparse matrix)
%
%   Example:
%

%   Copyright 2007-2007 Holger Heumann
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
  if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag; 
  else flags = zeros(nElements,1); end
  
  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices of current element
     %vertices
   
   vid=Mesh.Elements(i,:);
   Vertices=Mesh.Coordinates(vid,:);

   % Compute element contributions corresponding to the edges
    
   [Aloc1 Aloc2 Aloc3] = WEIGHT_WEDGE(Vertices);
  
   % Extract global edge numbers
    
   eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
           Mesh.Vert2Edge(vid(3),vid(1)) ...
           Mesh.Vert2Edge(vid(1),vid(2))];
       
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end
    
    Peori = diag([p1 p2 p3]); % scaling matrix taking into account orientations
    Aloc1 = Peori*Aloc1*Peori;
    Aloc2 = Peori*Aloc2*Peori;
    Aloc3 = Peori*Aloc3*Peori;
    
    Aloc1=Sigma(vid(1))*Aloc1+Sigma(vid(2))*Aloc2+Sigma(vid(3))*Aloc3;      
    
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(eidx,3);
    A(loc) = Aloc1(:);
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