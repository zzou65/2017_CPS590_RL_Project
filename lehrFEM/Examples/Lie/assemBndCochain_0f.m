function L = assemBndCochain_0f(Mesh, Vflags, FHandle)
% assemBndCochain_0f assemble Cochain representaion of zero forms on the
% boundary
%
%   assemBndCochain_0f calculates the cochain values for those vertices,
%   for which the sum of the Bdflags of the adjacent edges is equal to an
%   element in Vflags.
%
%   L = ASSEMBndCochain_0f(MESH, BdFlags, FHandel) .... 
%   Mesh:
%   BDFlags:    
%   FHandle:  vector valued function handle  
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemBndCochain_0f(Mesh,-1, Fhandle);
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%
% Initialize constants
% Assign output arguments
 
 Mesh=add_VBdFlags(Mesh);
 h=get_MeshWidth(Mesh);
 nC=size(Mesh.Coordinates,1);
 L=zeros(nC,1);
 for i=1:nC
     for j=Vflags
      if (Mesh.VBdFlags(i)==j)
       L(i)=FHandle(Mesh.Coordinates(i,:));
      end
     end
 end

return
   
  