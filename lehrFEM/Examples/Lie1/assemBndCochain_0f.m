function L = assemBndCochain_0f(Mesh, flags, FHandle)
% assemCochain_0f assemble Cochain of zero forms
%
%   A = ASSEMCochain_0f(MESH,FHandel) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMCochain_0f
%    D .... assembles the global matrix 
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassTwoD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 Mesh=add_VBdFlags(Mesh);
 h=get_MeshWidth(Mesh);
 nC=size(Mesh.Coordinates,1);
 L=zeros(nC,1);
 for i=1:nC
     for j=flags
      if (Mesh.VBdFlags(i)==j)
       L(i)=FHandle(Mesh.Coordinates(i,:));
      end
     end
 end

return
   
  