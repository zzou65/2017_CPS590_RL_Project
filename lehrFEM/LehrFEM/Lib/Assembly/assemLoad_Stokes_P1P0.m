function L = assemLoad_Stokes_P1P0(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_STOKES_P1P0 Assemble load vector for P1 elements.
%
%   L = ASSEMLOAD_STOKES_P1P0(MESH,QUADRULE,FHANDLE) assembles the global
%   load vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying the additional element
%                 information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_STOKES_P1P0(MESH,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function
%   handle FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)[x(:,1) x(:,2)];
%   L = assemLoad_Stokes_P1P0(Mesh,P7O6(),FHandle);
%
%   See also shap_LFE, assemLoad_Stokes_P2P0.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(2*nCoordinates+nElements+1,1);
  
  % Precompute shape functions
  
  N = shap_LFE(QuadRule.x);
  
  % Assemble element contributions

  for i = 1:nElements
  
    % Extract vertex and edge numbers
    
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
    
    offset = nCoordinates;
    L(vidx(1)+offset) = L(vidx(1)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,1))*det_BK;
    L(vidx(2)+offset) = L(vidx(2)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,2))*det_BK;
    L(vidx(3)+offset) = L(vidx(3)+offset) + ...
                        sum(QuadRule.w.*FVal(:,2).*N(:,3))*det_BK;
    
  end

return