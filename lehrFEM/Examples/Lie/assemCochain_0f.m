function L = assemCochain_0f(Mesh, FHandle)
% assemCochain_0f assemble Cochain of zero forms
%
%   A = ASSEMCochain_0f(MESH,FHandel) 
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemcochain_0f(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
L=FHandle(Mesh.Coordinates);

return
   
  