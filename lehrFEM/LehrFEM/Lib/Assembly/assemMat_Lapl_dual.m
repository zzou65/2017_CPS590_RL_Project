function varargout = assemMat_Lapl_dual(Mesh,EHandle,varargin)
% ASSEMMAT_LAPL_dual Assemble dual mixed FEM.
%
%   A = ASSEMMAT_LAPL_dual(MESH,EHANDLE) assembles the global matrix 
%   from the local element contributions given by the function handle 
%   EHANDLE and returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMAT_Lapl_dual(MESH,EHANDLE) assembles the global
%   matrix from the local element contributions given by the function 
%   handle EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%      A = assemMat_Lapl_dual(NewMesh,@STIMA_Lapl_dual,P7O6());

%   Copyright 2005-2006 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initilaize constants
  
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  I = zeros(16*nElements,1);
  J = zeros(16*nElements,1);
  A = zeros(16*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:16;
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to stiffness matrix
    
    eidx = [Mesh.Vert2Edge(vidx(2),vidx(3)) ...
            Mesh.Vert2Edge(vidx(3),vidx(1)) ...
            Mesh.Vert2Edge(vidx(1),vidx(2))];
    idx = [eidx ...
           nEdges+i ...
           ];
       
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vidx(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vidx(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vidx(1)),  p3 = 1;  else    p3 = -1;  end
    
   Peori = diag([p1 p2 p3 1]) ; % scaling matrix taking into account orientations
    Aloc = Peori*Aloc*Peori;        
    
           
    I(loc) = set_Rows(idx,4);
    J(loc) = set_Cols(idx,4);
    A(loc) = Aloc(:);
    loc = loc+16;
      
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