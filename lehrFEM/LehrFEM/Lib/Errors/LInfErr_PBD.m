function err = LInfErr_LFE(Mesh,u,FHandle,varargin)
% LINFERR_LFE Discretization error in Linf norm for linear finite elements
%             with parabolic boundary approximation.
%
%   ERR = LINFERR_PBD(MESH,U,FHANDLE) computes the discretization error
%   between the exact solution given by the function handle FHANDLE and the
%   finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of
%                 elements.
%    DELTA        P-by-1 matrix specifying the boundary correction term on
%                 every edge.
%
%   ERR = LINFERR_PBD(MESH,U,FHANDLE,FPARAM) also handles the variable
%   length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_PBD(Mesh,u,FHandle);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  Rot = [0 -1; 1 0];
  
  % Extract boundary edges, vertices and elements
  
  BndEdges = get_BdEdges(Mesh);
  BndNodes = unique([Mesh.Edges(BndEdges,1); Mesh.Edges(BndEdges,2)]);
  [rows,cols,vals] = find(Mesh.Edge2Elem(BndEdges,:));
  vals = sortrows([rows cols vals],1);
  BndElem = vals(:,3);
  
  % Compute mid points of straight edges
  
  IntEdges = setdiff(1:nEdges,BndEdges);
  XMid_s = (Mesh.Coordinates(Mesh.Edges(IntEdges,1),:) + ...
            Mesh.Coordinates(Mesh.Edges(IntEdges,2),:))/2;
  
  % Compute mid points of curved edges
   
  normal = Mesh.Coordinates(Mesh.Elements(BndElem,3),:) - ...
           Mesh.Coordinates(Mesh.Elements(BndElem,2),:);         
  normal = (normal*Rot)./(sqrt(sum(normal.^2,2))*ones(1,2));
  XMid_c = (Mesh.Coordinates(Mesh.Elements(BndElem,2),:) + ...
            Mesh.Coordinates(Mesh.Elements(BndElem,3),:))/2 + ...
           (Mesh.Delta(BndEdges)*ones(1,2)).*normal;
       
  % Measure error with respect to LInf norm   
     
  u_Vertices = u(1:nCoordinates);
  u_BndEdges = u(BndEdges+nCoordinates);
  u_IntEdges = u(IntEdges+nCoordinates);
  err = max([max(abs(u_Vertices-FHandle(Mesh.Coordinates,varargin{:}))) ...
             max(abs(u_BndEdges-FHandle(XMid_c,varargin{:}))) ...
             max(abs(u_IntEdges-FHandle(XMid_s,varargin{:})))]);

return