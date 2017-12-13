function [U,FreeDofs] = assemDir_W1F(Mesh,BdFlags,FHandle,QuadRule_1D,varargin)
% ASSEMDIR_W1F Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_W1F(MESH,BDFLAG,FHANDLE,QUADRULE_1D) incoporates 
%   the dirichlet boundary conditions with the data given by FHANDLE into 
%   the finite element solution U. The boundary condition is only enforced 
%   at the edges whose BdFlag is equal to the integer BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_W1F(MESH,BDFLAG,FHANDLE,QUADRULE_1D,FPARAM) 
%   also handles the variable length argument list FPARAM to the boundary 
%   data function FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each 
%                 boundary edge in the mesh.
%
%   FREEDOFS is a M-by-1 matrix specifying the vertices and edges with no
%   prescribed Dirichlet boundary data.  
%
%   Example:
%
%   [U,FreeDofs] = assemDir_W1F(Mesh,BdFlags,FHandle,gauleg(0,1,10));
%
%   See also get_BdEdges.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nEdges = size(Mesh.Edges,1);
  tmp = [];
  U = zeros(nEdges,1);
  nGuass = size(QuadRule_1D.w,1);
  
  for j = BdFlags
  
    % Extract Dirichlet nodes
  
    Loc = get_BdEdges(Mesh);
    DEdges = Loc(Mesh.BdFlags(Loc) == j);
    
    for i = 1:size(DEdges,1)

        P1 = Mesh.Coordinates(Mesh.Edges(DEdges(i),1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(DEdges(i),2),:);

        % Compute midpoints of all edges

        x = ones(nGuass,1)*P1+QuadRule_1D.x*(P2-P1);
        dS = ones(nGuass,1)*(P2-P1);
        Fval = FHandle(x,j,varargin{:});

        % Compute Dirichlet boundary conditions

        U(DEdges(i)) = sum(QuadRule_1D.w.*sum(Fval.*dS,2));

    end
                   
    % Collect Dirichlet nodes in temporary container
    
    tmp = [tmp ; DEdges];
    
  end
  
  % Compute set of free dofs
  
  FreeDofs = setdiff(1:nEdges,tmp);
  
return