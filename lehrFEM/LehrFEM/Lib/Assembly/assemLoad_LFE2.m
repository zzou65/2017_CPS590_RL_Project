function L = assemLoad_LFE2(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_LFE Assemble nodal FE contributions.
%
%   L = ASSEMLOAD_LFE2(MESH,QUADRULE,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle FHANDLE.
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
%   L = ASSEMLOAD_LFE2(COORDINATES,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemLoad_LFE2(Mesh,P7O6(),FHandle);
%
%   See also shap_LFE2.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(2*nCoordinates,1);
  
  % Precompute shape functions
  
  N = shap_LFE2(QuadRule.x);
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    
    % Compute load data
    
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal(:,1).*N(:,1))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal(:,1).*N(:,5))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal(:,1).*N(:,9))*det_BK;
      
    L(vidx(1)+nCoordinates) = L(vidx(1)+nCoordinates) + ...
                              sum(QuadRule.w.*FVal(:,2).*N(:,4))*det_BK; 
    L(vidx(2)+nCoordinates) = L(vidx(2)+nCoordinates) + ...
                              sum(QuadRule.w.*FVal(:,2).*N(:,8))*det_BK; 
    L(vidx(3)+nCoordinates) = L(vidx(3)+nCoordinates) + ...
                              sum(QuadRule.w.*FVal(:,2).*N(:,12))*det_BK;
                          
  end
  
return