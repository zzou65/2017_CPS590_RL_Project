function Mesh = refine_hp(Mesh,CNodes)
% REFINE_HP Adaptive refinement for hp-FEM.
%
%   MESH = REFINE_HP(MESH,CNODES) adaptively refines a mesh towards the
%   corner nodes CNODES using the longest edge bisection algorithm.
%   
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh. 
%    NEIGH       M-by-3 matrix specifying all elements sharing a common
%                edge with the current element. If an edge is part of the
%                boundary then NEIGH contains the boundary flag of that
%                edge.
%    OPP_VERT    M-by-3 matrix specifying the local number of the opposite
%                vertex of the element sharing an edge with the current
%                element.
%
%   Example:
%
%   Mesh = refine_hp(Mesh,CNodes);
%
%   See also refine_NVB.

%   Copyright 2006-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCNodes = size(CNodes,1);  % Number of corner nodes of the mesh 
  
  % Check elements for corner vertices
  
  MarkedElem = [];
  for i = 1:nCNodes
    for j = 1:3
      Idx = find(Mesh.Elements(:,j) == CNodes(i));
      MarkedElem = [MarkedElem; Idx]; 
    end
  end
  MarkedElem = unique(MarkedElem);
  
  % Do longest edge bisection with all marked elements
  
  Mesh = refine_LEB(Mesh,MarkedElem);
  
return