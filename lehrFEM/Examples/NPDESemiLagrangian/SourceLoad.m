function L = SourceLoad(Mesh,QuadRule,FHandle,t)
% SourceLoad Computes a Source Function Vector for a Mesh at Time instant t 
% using function handle Fhandle and quadrature rule QuadRule.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%
%   By Rajdeep Deb
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
    
    FVal = FHandle(x(:,1),x(:,2),t);
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1))*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3))*det_BK;
      
  end
  
return