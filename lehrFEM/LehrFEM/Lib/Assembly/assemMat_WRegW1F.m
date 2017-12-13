function varargout = assemMat_WRegW1F(Mesh,EHandle,varargin)
% ASSEMMAT_WREGW1F Assemble WREG W1F FE contributions.
%
%   A = ASSEMMAT_WREGW1F(MESH,EHANDLE) assembles the global matrix from the
%   local element contributions given by the function handle EHANDLE and
%   returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_WREGW1F(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   [I,J,A] = ASSEMMAT_WREGW1F(MESH,EHANDLE) assembles the global matrix 
%   from the local element contributions given by the function handle 
%   EHANDLE and returns the matrix in an array representation.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   EHandle = @STIMA_WReg_W1F;
%   A = assemMat_WRegW1F(Mesh,EHandle);
%  
%   See also SET_ROWS, SET_COLS.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  
  % Preallocate memory
  
  I = zeros(9*nElements,1);
  J = zeros(9*nElements,1);
  A = zeros(9*nElements,1);
  
  % Assemble element contributions
  
  loc = 1:9;
  
  for i = 1:nElements
    
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
      
    % Compute element contributions
    
    Aloc = EHandle(Vertices,Mesh.ElemFlag(i),P7O6(),varargin{:});
   
    
    % Extract global edge numbers
    
    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
       
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vidx(2))
        p1 = 1;
    else
        p1 = -1;
    end
    
    if(Mesh.Edges(eidx(2),1)==vidx(3))
        p2 = 1;
    else
        p2 = -1;
    end
    
    if(Mesh.Edges(eidx(3),1)==vidx(1))
        p3 = 1;
    else
        p3 = -1;
    end
    
    Peori = diag([p1 p2 p3]);
    Aloc = Peori*Aloc;
       
    % Add contributions to stiffness matrix
    
    I(loc) = set_Rows(eidx,3);
    J(loc) = set_Cols(vidx,3);
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