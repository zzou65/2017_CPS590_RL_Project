function L = assemCochain_1f(Mesh, FHandle, QuadRule_1D,varargin)
% assemCochain_1f assemble Cochain of one forms
%
%   L = ASSEMCochain_1f(MESH,FHandel) .... 
%
%
%   Copyright 2007-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
 nEdges=size(Mesh.Edges,1);
 L = zeros(nEdges,1);
 nGuass = size(QuadRule_1D.w,1);
 
 for i = 1:nEdges

        P1 = Mesh.Coordinates(Mesh.Edges(i,1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(i,2),:);

        % Compute midpoints of all edges

        x = ones(nGuass,1)*P1+QuadRule_1D.x*(P2-P1);
        dS = ones(nGuass,1)*(P2-P1);
        Fval = FHandle(x,0,varargin{:});

        % Compute Dirichlet boundary conditions

        L(i) = sum(QuadRule_1D.w.*sum(Fval.*dS,2));

 end
  
return
   
  