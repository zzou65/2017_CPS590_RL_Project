function L = assemRotContr_Bnd_W1F(Mesh, BdFlags, FHandle, VHandle, QuadRule_1D,varargin)
% assemRotContr_Bnd_W1F assembles Neumann-like boundary data for weak
% treatment of curl (v x u) terms
%
%   Copyright 2007-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
 nEdges=size(Mesh.Edges,1);
 
 L = zeros(nEdges,1);
 nGauss = size(QuadRule_1D.w,1);
 
 for i = 1:nEdges
     for j=BdFlags
      if (Mesh.BdFlags(i)==j)
        P1 = Mesh.Coordinates(Mesh.Edges(i,1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(i,2),:);

        % Compute quadrature points
        x  = ones(nGauss,1)*P1+QuadRule_1D.x*(P2-P1);
        dS = norm(P2-P1);
        Uval = FHandle(x,0,varargin{:});
        Vval = VHandle(x);
        Uval = [Uval(:,2), -Uval(:,1)]; 

        % conditions
        L(i) = sum(QuadRule_1D.w.*sum(Vval.*Uval,2))*dS;
        
      end
     end
 end
  
return
   
  