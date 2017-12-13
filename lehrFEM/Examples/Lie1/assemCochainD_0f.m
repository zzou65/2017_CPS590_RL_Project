function L = assemCochainD_0f(Mesh, FHandle)
% assemCochainD_0f assemble Cochain of zero forms on dual grid
%
%   A = ASSEMCochainD_0f(MESH,FHandel) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMCochainD_0f
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
 
  b_Mesh=refine_BAR(Mesh);
  b_Mesh = add_Edge2Elem(b_Mesh);
  b_Mesh = add_Patches(b_Mesh);
  b_nElements=size(b_Mesh,1);
  b_U=zeros(b_nElements,1);
  nCoordinates=size(Mesh.Coordinates);
  L=zeros(nCoordinates,1);
  b_L=assemCochain0f(b_Mesh,FHandle);
  for i=1:nCoordinates
      b_AdjElements=b_Mesh.AdjElements(i,:);
      b_AdjElements=setdiff(b_AdjElements,0);
      b_nAdjElements=size(b_AdjElements,2);
      L=sum
  end
L=FHandle(Mesh.Coordinates);

return
   
  