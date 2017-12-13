function L = assemLoad_Lapl_dual(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_Lapl_dual assemble load vector for dual laplace.
%
%   L = ASSEMLOAD_LOAD_Lapl_dual(MESH,QUADRULE,FHANDLE) assembles the global
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
%    X Abscissae of the Gauss quadrature
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1) x(:,2)];
%   L = assemLoad_Lapl_dual(Mesh,P7O6(),FHandle);
%
%   See also shap_CR.

%   Copyright 2005-2006 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(nEdges+nElements,1);
  
  
  % Assemble element contributions
 
  eidx = zeros(1,3);
  for i = 1:nElements
  
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = -FHandle(x,Mesh.ElemFlag(i),varargin{:});

    % Add contributions to global load vector
    
    L(nEdges+i) = L(nEdges+i) + sum(QuadRule.w.*FVal.*2);
    
  end

return