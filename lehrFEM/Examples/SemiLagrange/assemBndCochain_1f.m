function L = assemBndCochain_1f(Mesh, BdFlags, FHandle, QuadRule_1D)
% assemCochain_1f assemble Cochain of one forms
%
%   L = ASSEMCochain_1f(MESH,FHandel) .... 
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
 
 nEdges=size(Mesh.Edges,1);
 
 L = zeros(nEdges,1);
 nGuass = size(QuadRule_1D.w,1);
 
 for i = 1:nEdges
     for j=BdFlags
      if (Mesh.BdFlags(i)==j)
        P1 = Mesh.Coordinates(Mesh.Edges(i,1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(i,2),:);

        % Compute midpoints of all edges

        x  = ones(nGuass,1)*P1+QuadRule_1D.x*(P2-P1);
        dS = ones(nGuass,1)*(P2-P1);
        Fval = FHandle(x);

        % Compute Dirichlet boundary conditions

        L(i) = sum(QuadRule_1D.w.*sum(Fval.*dS,2));
      end
     end
 end
  
return
   
  