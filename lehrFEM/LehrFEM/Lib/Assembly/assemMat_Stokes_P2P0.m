function varargout = assemMat_Stokes_P2P0(Mesh,EHandle,varargin)
% ASSEMMAT_STOKES_P2P0 Assemble Galerking matrix for P2-P0 finite element discretization
% of Stokes problem
%
%   A = ASSEMMAT_STOKES_P2P0(MESH,EHANDLE) assembles the global matrix 
%   from the local element contributions given by the function handle 
%   EHANDLE and returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_STOKES_P2P0(MESH,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during 
%   the assembly process. 
%
%   [I,J,A] = ASSEMMAT_STOKES_P2P0(MESH,EHANDLE) assembles the global 
%   matrix from the local element contributions given by the function 
%   handle EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or matrix specifying the elements of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   EHandle = @STIMA_Stokes_P2P0;
%   A = assemMat_Stokes_P2P0(Mesh,@STIMA_Stokes_P2P0,NU,P7O6());

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initilaize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  I = zeros(196*nElements,1);
  J = zeros(196*nElements,1);
  A = zeros(196*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:196;
  for i = 1:nElements
    
    % Extract vertices
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute 14x14 element matrix
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to stiffness matrix
    eidx = [Mesh.Vert2Edge(vidx(1),vidx(2)) ...
            Mesh.Vert2Edge(vidx(2),vidx(3)) ...
            Mesh.Vert2Edge(vidx(3),vidx(1))];
    idx = [vidx ...
           eidx+nCoordinates ...
           vidx+nCoordinates+nEdges ...
           eidx+2*nCoordinates+nEdges ...
           i+2*(nEdges+nCoordinates) ...
           2*(nEdges+nCoordinates)+nElements+1];
    I(loc) = set_Rows(idx,14);
    J(loc) = set_Cols(idx,14);
    A(loc) = Aloc(:);
    loc = loc+196;
  end
  
  % Assign output arguments for creation of sparse matrix
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);  
  end
  
return