function L = assemLoad_W1F2nd(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_W1F Assemble W1F FE contributions.
%
%   L = ASSEMLOAD_W1F(MESH,QUADRULE,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle EHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_W1F(COORDINATES,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1).^2 x(:,2).^2];
%   L = assemLoad_W1F(Mesh,P7O6(),FHandle);
%
%   See also shap_W1F.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Check for element flags
  if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag; 
  else flags = zeros(nElements,1); end
  
  % Preallocate memory
  
  L = zeros(2*nEdges,1);
  
  % Precompute shape functions
  
  N = shap_W1F2nd(QuadRule.x);
  
  % Assemble element contributions
  
  eidx = zeros(1,3);
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    eidx(1) = Mesh.Vert2Edge(vidx(2),vidx(3));
    eidx(2) = Mesh.Vert2Edge(vidx(3),vidx(1));
    eidx(3) = Mesh.Vert2Edge(vidx(1),vidx(2));
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    TK = transpose(inv(BK));
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = FHandle(x,flags(i),varargin{:});
    
    % Determine the orientation
    
    if(Mesh.Edges(eidx(1),1)==vidx(2))
        p1 = 1;
    else
        p1 = -1;
    end
    
    if(Mesh.Edges(eidx(2),1)==vidx(3))
        p2 = 1;
    else
        p2 = -1;
    end
    
    if(Mesh.Edges(eidx(3),1)==vidx(1))
        p3 = 1;
    else
        p3 = -1;
    end
    
    % Add contributions to global load vector
        
    L(eidx(1)) = L(eidx(1)) + sum(QuadRule.w.*sum(FVal.*([N(:,1) N(:,2)]*TK),2))*det_BK*p1;
    L(eidx(2)) = L(eidx(2)) + sum(QuadRule.w.*sum(FVal.*([N(:,3) N(:,4)]*TK),2))*det_BK*p2;
    L(eidx(3)) = L(eidx(3)) + sum(QuadRule.w.*sum(FVal.*([N(:,5) N(:,6)]*TK),2))*det_BK*p3;
    L(eidx(1)+nEdges) = L(eidx(1)+nEdges) + sum(QuadRule.w.*sum(FVal.*([N(:,7) N(:,8)]*TK),2))*det_BK;
    L(eidx(2)+nEdges) = L(eidx(2)+nEdges) + sum(QuadRule.w.*sum(FVal.*([N(:,9) N(:,10)]*TK),2))*det_BK;
    L(eidx(3)+nEdges) = L(eidx(3)+nEdges) + sum(QuadRule.w.*sum(FVal.*([N(:,11) N(:,12)]*TK),2))*det_BK;
      
  end
  
return