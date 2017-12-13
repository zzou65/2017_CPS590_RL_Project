
function L = assemCochainD_2f(Mesh, FHandle,QuadRule)
% assemCochainD_2f assemble Cochain of two forms on dual grid
%
%   A = ASSEMCochainD_2f(MESH,FHandel) .... and
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

% barycentric refinement 

  b_Mesh=refine_BAR(Mesh);
  b_Mesh = add_Edge2Elem(b_Mesh);
  b_Mesh = add_Patches(b_Mesh);
  b_nElements=size(b_Mesh,1);
  b_L=zeros(b_nElements,1);
  
% cochain representation on finer mesh  
 
  b_L=assemCochain_2f(b_Mesh,FHandle,QuadRule);
  
% Assign output arguments
  
  nCoordinates=size(Mesh.Coordinates,1);
  L=zeros(nCoordinates,1);
  
% Distribution to duals of Meshvertices (cells)  
  for i=1:nCoordinates
      b_AdjElements=b_Mesh.AdjElements(i,:);
      b_AdjElements=setdiff(b_AdjElements,0);
      L(i)=sum(b_L(b_AdjElements));
  end
return
   
  