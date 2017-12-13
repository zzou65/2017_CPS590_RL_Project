function [L] = SemiLagConstantVelocity_LFE(Mesh,U0,tau,Convectionhandle,Dhandle,StoreElement,QuadRule)
% SemiLagConstantVelocity_LFE Asseble Load Vector Corresponding to previous time step 
% and lagrangian tracking position with stationary convection speed.
%
%   L = SemiLagrangian_LFE(Mesh,U0,tau,Convectionhandle,Dhandle) assembles  
%   Load Vector based on Convectionhandle
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   by Rajdeep Deb
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  
  %Load Vector U(x-tauV)*V(x) are integrated based on P3O2() Gaussian
  %Quadrature
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    flags = Mesh.ElemFlag; 
  else
    flags = zeros(nElements,1);
  end
  
     Mesh.ElemFlag = flags;

  N = shap_LFE(QuadRule.x); 
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
    
    %construct load corresponding to previous time step
    
    %construct LFE using U and linear functions
       
    [WU] = computeLoadIntegral_ConstantVelocity(U0,QuadRule,Vertices,Mesh,tau,Convectionhandle,Dhandle,StoreElement(i,:));
    
    bK = Mesh.Coordinates(idx(1),:);
    BK = [Mesh.Coordinates(idx(2),:)-bK; Mesh.Coordinates(idx(3),:)-bK];
    det_BK = abs(det(BK));
      
    
    L(idx(1)) = L(idx(1)) + sum(WU.*N(:,1))*det_BK;
    L(idx(2)) = L(idx(2)) + sum(WU.*N(:,2))*det_BK;
    L(idx(3)) = L(idx(3)) + sum(WU.*N(:,3))*det_BK;
    
  end 