function L = assemLoad_BFE(Mesh,QuadRule,FHandle,varargin)
% ASSEMLOAD_BFE Assemble bilinear FE contributions.
%
%   L = ASSEMLOAD_BFE(MESH,QUADRULE,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle EHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-4 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_BFE(COORDINATES,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemLoad_BFE(Mesh,TProd(gauleg(0,1,2)),FHandle);
%
%   See also shap_BFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initilaize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertex numbers
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element mapping
    
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    P4 = Mesh.Coordinates(vidx(4),:);
  
    N = shap_BFE(QuadRule.x);
    grad_N = grad_shap_BFE(QuadRule.x);
  
    z1 = P1(1)*grad_N(:,1:2)+P2(1)*grad_N(:,3:4) + ...
         P3(1)*grad_N(:,5:6)+P4(1)*grad_N(:,7:8);
    z2 = P1(2)*grad_N(:,1:2)+P2(2)*grad_N(:,3:4) + ...
         P3(2)*grad_N(:,5:6)+P4(2)*grad_N(:,7:8);
    det_DPhi_K = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));
  
    x = N(:,1)*P1+N(:,2)*P2+N(:,3)*P3+N(:,4)*P4;
  
    % Compute load data
  
    FVal = FHandle(x,Mesh.ElemFlag(i),varargin{:});
    
    % Add contributions to global load data
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal.*N(:,1).*det_DPhi_K); 
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal.*N(:,2).*det_DPhi_K);
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal.*N(:,3).*det_DPhi_K); 
    L(vidx(4)) = L(vidx(4)) + sum(QuadRule.w.*FVal.*N(:,4).*det_DPhi_K);
      
  end
  
return