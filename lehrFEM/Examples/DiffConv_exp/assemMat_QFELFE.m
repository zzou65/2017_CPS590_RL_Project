function varargout = assemMat_QFELFE(Mesh, F_Handle, varargin)
% assemMat_ContrOne Assemble topological gradient.
%
%   A = assemMat_Contr1f(Mesh, vHandle, varargin) 
%   A = ASSEMMat_Contr1f(MESH,  vHandle) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Contr1f(MESH, vHandle) .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)]
%   A = assemMat_Contr1f(Mesh,V_Handle);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges =size(Mesh.Edges,1);
  nElements=size(Mesh.Elements,1);
   
  % Preallocate memory
  
  I = zeros(18*nElements,1);
  J = zeros(18*nElements,1);
  A = zeros(18*nElements,1);
  
  % local MASS-Matrix, exact
%   BS=1/120*[2 -1 -1 8 4 8; ...
%                             -1 2 -1 8 8 4; ...
%                             -1 -1 2 4 8 8];
  
  % look for upwind nodes on each triangle
  loc=1:18;
  for i = 1:nElements
     B=zeros(3,6);   
     % Vertices
     vid = Mesh.Elements(i,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
          
     % Extract global edge numbers
     eidx = [Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1))];
     idx=eidx+nCoordinates;   
        
     % Compute element mapping
     
     bK = a1;
     BK = [a2-bK; ...
        a3-bK];
     det_BK = abs(det(BK));
     inv_BK = inv(BK); 
     
     % Element Matrices

     B=F_Handle([a1;a2;a3,],0,varargin{:}); 
     
     I(loc) = set_Rows(vid ,6);
     J(loc) = set_Cols([vid idx],3);
     A(loc) = B(:);
     loc = loc+18;
  end
  
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nCoordinates,nCoordinates+nEdges);      
  end
  
return