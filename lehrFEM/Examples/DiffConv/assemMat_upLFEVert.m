function varargout = assemMat_upLFEVert(Mesh, VHandle,weights)
% assemMat_UpLFE assembles Matrices M,UP and C to build upwind scheme 
%    M: diagonal Massmatrix
%    Up: Upwind part, using DOFS on element boundary
%    C:   Central part, using quadrature points inside element
%
%   A = assemMat_UpLFE(Mesh,Vhandle) 
%   A = ASSEMMat_UpLFE(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_UpLFE .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%    
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassUpQFE(Mesh,VHandle);
%  
%   Copyright 2008-2008 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
  nF=size(Mesh.Elements,1);
  nE=size(Mesh.Edges,1);
  nC=size(Mesh.Coordinates,1);
  
  % Preallocate memory
  
  % upwinding part
  I_UP = zeros(9*nF,1);
  J_UP = zeros(9*nF,1);
  UP = zeros(9*nF,1);
    
  MD=zeros(nE,1);
  
  for f=1:nF
     % Vertices
     vid = Mesh.Elements(f,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % Element Mass
     BK = [a2-a1;a3-a1];
     det_BK = abs(det(BK));
     
     MD(vid)=MD(vid)+det_BK/6*[1 ;1; 1];
  end

  loc = 1:9;
  
  for (f = 1:nF)
       
     % Vertices
     vid = Mesh.Elements(f,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
 
     % Compute element contributions
     UP_loc=STIMA_UPLFEVert([a1; a2; a3],VHandle, MD(vid));
     
     % Add contributions to global matrix
     I_UP(loc) = set_Rows(vid,3);
     J_UP(loc) = set_Cols(vid,3);
     UP(loc) = UP_loc(:);
        
     loc = loc+9;
        
  end
  
%  if(nargout > 1)
    varargout{1} = sparse(I_UP,J_UP,UP) ;
%   else
%     varargout{1} = [sparse(I,I,MD);      
%   end  
   
  return