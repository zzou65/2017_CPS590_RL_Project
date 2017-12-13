function L = assemCochain_0f(Mesh, FHandle)
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
 
L=FHandle(Mesh.Coordinates);

return
   
  