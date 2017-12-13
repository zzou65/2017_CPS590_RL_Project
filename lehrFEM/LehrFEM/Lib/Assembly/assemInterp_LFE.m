function L = assemInterp_LFE(Mesh,QuadRule,FHandle,GDHandle,varargin)
% ASSEMINTERP_LFE Assemble linear FE contributions.
%
%   L = ASSEMINTERP_LFE(MESH,QUADRULE,FHANDLE,GDHANDLE) assembles the global load 
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
%   L = ASSEMINTERP_LFE(COORDINATES,QUADRULE,FHANDLE,GDHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemInterp_LFE(Mesh,P7O6(),FHandle,GDHandle);
%
%   See also shap_LFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Precompute shape functions
  
  N = shap_LFE(QuadRule.x);
  GN_R = grad_shap_LFE(QuadRule.x);
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    
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
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1))*det_BK;
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*sum(GDval.*GN(:,1:2),2))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*sum(GDval.*GN(:,3:4),2))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*sum(GDval.*GN(:,5:6),2))*det_BK;
      
  end
  
return