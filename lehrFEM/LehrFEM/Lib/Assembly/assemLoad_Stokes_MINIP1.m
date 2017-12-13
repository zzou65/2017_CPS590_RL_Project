function L = assemLoad_Stokes_MINIP1(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_STOKES_MINIP1 Assemble load vector for MINI elements.
%
%   L = ASSEMLOAD_STOKES_MINIP1(MESH,QUADRULE,FHANDLE) assembles the global
%   load vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by_3 or N-by-4 matrix specifying the elements of the
%                 mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_STOKES_MINIP1(MESH,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function
%   handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1) x(:,2)];
%   L = assemLoad_Stokes_MINIP1(Mesh,P7O6(),FHandle);
%   
%   See also shap_MINI.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(3*nCoordinates+2*nElements+1,1);
  
  % Precompute shape functions
  
  N = shap_MINI(QuadRule.x); 
  
  % Assemble element contributions
  
  eidx = zeros(1,3);
  for i = 1:nElements
  
    % Extract vertex numbers
    
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
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal(:,1).*N(:,2))*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal(:,1).*N(:,3))*det_BK;
    L(i+nCoordinates) = L(i+nCoordinates) + ...
                        sum(QuadRule.w.*FVal(:,1).*N(:,4))*det_BK;
    
    offset = nCoordinates+nElements;
    L(vidx(1)+offset) = L(vidx(1)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,1))*det_BK;
    L(vidx(2)+offset) = L(vidx(2)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,2))*det_BK;
    L(vidx(3)+offset) = L(vidx(3)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,3))*det_BK;
    L(i+offset+nCoordinates) = L(i+offset+nCoordinates) + ...
                               sum(QuadRule.w.*FVal(:,2).*N(:,4))*det_BK;
      
  end

return