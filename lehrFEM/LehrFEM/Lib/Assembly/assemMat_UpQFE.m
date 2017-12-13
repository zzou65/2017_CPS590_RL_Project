function varargout = assemMat_UpQFE(Mesh, VHandle,weights)
% assemMat_UpQFE assembles Matrices M,UP and C to build upwind scheme M*Up+C 
%    M: diagonal Massmatrix
%    Up: Upwind part, using DOFS on element boundary
%    C:   Central part, using quadrature points inside element
%
%   A = assemMat_UpQFE(Mesh,Vhandle) 
%   A = ASSEMMat_UpQFE(MESH) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_UpQFE .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%    
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassUpQFE(Mesh,VHandle);
%  
%   Copyright 2007-2008 Holger Heumann
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
  I_UP = zeros(36*nF,1);
  J_UP = zeros(36*nF,1);
  UP = zeros(36*nF,1);
  
  % part for barycentric quadrature points
  I_C = zeros(36*nF,1);
  J_C = zeros(36*nF,1);
  C = zeros(36*nF,1);
  
  MD=zeros(nC+nE,1);
  I_MD=1:(nC+nE);
  
  loc = 1:36;
  
  for (f = 1:nF)
       
     % Vertices
     vid = Mesh.Elements(f,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % Compute element contributions
     [MD_loc, UP_loc, C_loc]=STIMA_UPQFE([a1; a2; a3],VHandle,weights);
     
     % Extract global edge numbers
    
     idx = [vid ...
           Mesh.Vert2Edge(Mesh.Elements(f,1),Mesh.Elements(f,2))+nC ...
           Mesh.Vert2Edge(Mesh.Elements(f,2),Mesh.Elements(f,3))+nC ...
           Mesh.Vert2Edge(Mesh.Elements(f,3),Mesh.Elements(f,1))+nC];

     % Add contributions to global matrix
      
     MD(idx)=MD(idx)+MD_loc';
    
     I_UP(loc) = set_Rows(idx,6);
     J_UP(loc) = set_Cols(idx,6);
     UP(loc) = UP_loc(:);
     
     I_C(loc) = set_Rows(idx,6);
     J_C(loc) = set_Cols(idx,6);
     C(loc) = C_loc(:);
     
     loc = loc+36;
    
     
  end
  
%  if(nargout > 1)
    varargout{1} = sparse(I_MD,I_MD,MD);
    varargout{2} = sparse(I_UP,J_UP,UP) ;
    varargout{3} = sparse(I_C,J_C,C);
%   else
%     varargout{1} = [sparse(I,I,MD);      
%   end  
   
  return
