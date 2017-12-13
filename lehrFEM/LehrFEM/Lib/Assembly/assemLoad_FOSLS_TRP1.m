function L = assemLoad_FOSLS_TRP1(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_FOSLS_TRP1 Assemble load vector for Crouzeix-Raviart elemets.
%
%   L = ASSEMLOAD_FOSLS_TRP1(MESH,QUADRULE,FHANDLE) assembles the global
%   load vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying the additional element
%                 information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_FOSLS_TRP1(MESH,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1) x(:,2)];
%   L = assemLoad_FOSLS_TRP1(Mesh,P7O6(),FHandle);
%
%   See also shap_CR.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(nEdges+nCoordinates,1);
  
  
  % Assemble element contributions
  
  eidx = zeros(1,3);
  for i = 1:nElements
  
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(2),vidx(3));
    eidx(2) = Mesh.Vert2Edge(vidx(3),vidx(1));
    eidx(3) = Mesh.Vert2Edge(vidx(1),vidx(2));
      
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
    
    
    if(Mesh.Edges(eidx(1),1)==vidx(2)),  p1 = 1;  else    p1 = -1;  end
    if(Mesh.Edges(eidx(2),1)==vidx(3)),  p2 = 1;  else    p2 = -1;  end
    if(Mesh.Edges(eidx(3),1)==vidx(1)),  p3 = 1;  else    p3 = -1;  end
    
    % Add contributions to global load vector
    
    L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*FVal(:,1).*2.*(p1));
    L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*FVal(:,1).*2.*(p2));
    L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*FVal(:,1).*2.*(p3));
    
  end

return