function L = assemNeu_LFV(mesh,bdFlags,L,quadRule,fHandle,varargin)
% ASSEMNEU_LFV Neumann boundary conditions.
%
%   L = ASSEMNEU_LFV(MESH,BDFLAG,L,QUADRULE,FHANDLE) incoporates the
%   Neumann boundary conditions with the data given by FHANDLE into the
%   right hand side load vector L. The boundary condition is only enforced
%   at edges whose BdFlag is equal to the integer BDFLAG. The 1D struct
%   QUADRULE is used to do the numerical integration along the edges.
%
%   L = ASSEMNEU_LFV(MESH,BDFLAG,L,QUADRULE,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the boundary data function
%   FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first column
%                 specifies the left hand side element where the second column
%                 specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of elements. 
%
%   Example:
%
%   L = assemNeu_LFV(Mesh,BdFlags,L,QuadRule,FHandle);
%
%   See also GET_BDEDGES, SHAP_LFE.

%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    bdEdges = get_BdEdges(mesh);	

    for f=bdFlags
	edges = bdEdges(mesh.BdFlags(bdEdges) == f);
	for e=edges'
	    i = mesh.Edges(e,1);
	    j = mesh.Edges(e,2);
	    vi = mesh.Coordinates(i,:);
	    vj = mesh.Coordinates(j,:);
	    mp = mesh.MidPoints(e,:);

	    inti = 0; intj = 0;
	    for a=1:length(quadRule.w)
		    inti = inti + quadRule.w(a)*fHandle(vi+quadRule.x(a)*(mp-vi),varargin)*norm(mp-vi);
		    intj = intj + quadRule.w(a)*fHandle(vj+quadRule.x(a)*(mp-vj),varargin)*norm(mp-vj);
	    end
	    L(i) = L(i) + inti;
	    L(j) = L(j) + intj;
	end
    end

end
