function varargout = assemMat_Stokes_TH(Mesh,EHandle,varargin)
% ASSEMMAT_STOKES_TH Assemble Taylor-Hood FE contributions.
%
%   A = ASSEMMAT_STOKES_TH(MESH,EHANDLE) assembles the global matrix 
%   from the local element contributions given by the function handle 
%   EHANDLE and returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_STOKES_TH(MESH,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during 
%   the assembly process. 
%
%   [I,J,A] = ASSEMMAT_STOKES_TH(MESH,EHANDLE) assembles the global 
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
%   EHandle = @STIMA_Stokes_TH;
%   A = assemMat_Stokes_TH(Mesh,@STIMA_Stokes_TH,NU,P7O6());

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initilaize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  I = zeros(256*nElements,1);
  J = zeros(256*nElements,1);
  A = zeros(256*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:256;
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to stiffness matrix
    
    eidx = [Mesh.Vert2Edge(vidx(1),vidx(2)) ...
            Mesh.Vert2Edge(vidx(2),vidx(3)) ...
            Mesh.Vert2Edge(vidx(3),vidx(1))];
    idx = [vidx ...
           eidx+nCoordinates ...
           vidx+nCoordinates+nEdges ...
           eidx+2*nCoordinates+nEdges ...
           vidx+2*(nEdges+nCoordinates) ...
           2*nEdges+3*nCoordinates+1];
    I(loc) = set_Rows(idx,16);
    J(loc) = set_Cols(idx,16);
    A(loc) = Aloc(:);
    loc = loc+256;
      
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