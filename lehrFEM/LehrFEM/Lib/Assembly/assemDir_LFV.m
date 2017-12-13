function [U,FreeDofs] = assemDir_LFV(mesh,bdFlags,fHandle,varargin)
% ASSEMDIR_LFV Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_LFV(MESH,BDFLAG,FHANDLE) incoporates the
%   Dirichlet boundary conditions with the data given by FHANDLE into the
%   finite volume solution U. The boundary condition is only enforced at
%   the vertices of the edges whose BdFlag is equal to the integer BDFLAG.
%
%   [U,FREEDOFS] = ASSEMDIR_LFE(MESH,BDFLAG,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the boundary data
%   function FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each
%                 boundary edge in the mesh.
%
%   FREEDOFS is a M-by-1 matrix specifying the vertices with no prescribed
%   Dirichlet boundary data.  
%
%   Example:
%
%   [U,FreeDofs] = assemDir_LFV(Mesh,BdFlags,FHandle);
%
%   See also get_BdEdges.

%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

	% Initialize

	nCoords = size(mesh.Coordinates,1);
	tmp = [];
	U = zeros(nCoords,1);

	for j=bdFlags
		% Get nodes
		locs = get_BdEdges(mesh);
		dEdges = locs(mesh.BdFlags(locs) == j);
		dNodes = unique([mesh.Edges(dEdges,1); mesh.Edges(dEdges,2)]);
	
		% Compute U-vector
		for i=dNodes'
			U(i) = fHandle(mesh.Coordinates(i,:),j,varargin{:});
		end

		% Collect nodes
		tmp = [tmp; dNodes];
	end

	FreeDofs = setdiff(1:nCoords,tmp);

return
