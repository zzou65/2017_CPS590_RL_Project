function L = assemInterp_QFE(Mesh,QuadRule,FHandle,GDHandle,varargin)
% ASSEMINTERP_QFE Assemble quadratic FE contributions.
%
%   L = ASSEMINTERP_QFE(MESH,QUADRULE,FHANDLE,GDHANDLE) assembles the global load 
%   vector for the load data given by the function handle EHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
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
%   L = ASSEMINTERP_QFE(COORDINATES,QUADRULE,FHANDLE,GDHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemInterp_QFE(Mesh,P7O6(),FHandle,GDHandle);
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates+nEdges,1);

  % Precompute shape functions
  
  N = shap_QFE(QuadRule.x);
  GN_R = grad_shap_QFE(QuadRule.x);
  eidx = zeros(1,3);
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(1),vidx(2)) + nCoordinates;
    eidx(2) = Mesh.Vert2Edge(vidx(2),vidx(3)) + nCoordinates;
    eidx(3) = Mesh.Vert2Edge(vidx(3),vidx(1)) + nCoordinates;
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    inv_trans_BK = inv(transpose(BK));
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
    GDval = GDHandle(x,Mesh.ElemFlag(i),varargin{:});
    GN(:,1:2) = GN_R(:,1:2)*inv_trans_BK;
    GN(:,3:4) = GN_R(:,3:4)*inv_trans_BK;
    GN(:,5:6) = GN_R(:,5:6)*inv_trans_BK;
    GN(:,7:8) = GN_R(:,7:8)*inv_trans_BK;
    GN(:,9:10) = GN_R(:,9:10)*inv_trans_BK;
    GN(:,11:12) = GN_R(:,11:12)*inv_trans_BK;
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1))*det_BK;
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*sum(GDval.*GN(:,1:2),2))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*sum(GDval.*GN(:,3:4),2))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*sum(GDval.*GN(:,5:6),2))*det_BK;
    L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*FVal.*N(:,4))*det_BK;
    L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*sum(GDval.*GN(:,7:8),2))*det_BK;
    L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*FVal.*N(:,5))*det_BK;
    L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*sum(GDval.*GN(:,9:10),2))*det_BK;
    L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*FVal.*N(:,6))*det_BK;
    L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*sum(GDval.*GN(:,11:12),2))*det_BK;
      
  end
  
return