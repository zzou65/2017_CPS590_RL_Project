function L = assemLoad_Stokes_P2P0(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_STOKES_P2P0 Assemble load vector for P2 elements.
%
%   L = ASSEMLOAD_STOKES_P2P0(MESH,QUADRULE,FHANDLE) assembles the global
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
%   L = ASSEMLOAD_STOKES_P2P0(MESH,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function
%   handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1) x(:,2)];
%   L = assemLoad_Stokes_P2P0(Mesh,P7O6(),FHandle);
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(2*(nCoordinates+nEdges)+nElements+1,1);
  
  % Precompute shape functions
  
  N = shap_QFE(QuadRule.x);
  
  % Assemble element contributions
  
  eidx = zeros(1,3);
  for i = 1:nElements
  
    % Extract vertex and edge numbers
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(1),vidx(2)) + nCoordinates;
    eidx(2) = Mesh.Vert2Edge(vidx(2),vidx(3)) + nCoordinates;
    eidx(3) = Mesh.Vert2Edge(vidx(3),vidx(1)) + nCoordinates;
      
    % Compute element mapping 
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal(:,1).*N(:,1))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal(:,1).*N(:,2))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal(:,1).*N(:,3))*det_BK;
    L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*FVal(:,1).*N(:,4))*det_BK;
    L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*FVal(:,1).*N(:,5))*det_BK;
    L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*FVal(:,1).*N(:,6))*det_BK;
    
    offset = nCoordinates + nEdges;
    L(vidx(1)+offset) = L(vidx(1)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,1))*det_BK;
    L(vidx(2)+offset) = L(vidx(2)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,2))*det_BK;
    L(vidx(3)+offset) = L(vidx(3)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,3))*det_BK;
    L(eidx(1)+offset) = L(eidx(1)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,4))*det_BK;
    L(eidx(2)+offset) = L(eidx(2)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,5))*det_BK;
    L(eidx(3)+offset) = L(eidx(3)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,6))*det_BK;
    
  end

return