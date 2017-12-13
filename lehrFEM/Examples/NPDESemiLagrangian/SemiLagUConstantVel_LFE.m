function [L] = SemiLagUConstantVel_LFE(Mesh,U0,tau,Convectionhandle,Dhandle,StoreElement)
% SemiLagConstantVelocity_LFE Asseble Load Vector Corresponding to previous time step 
% and lagrangian tracking position with stationary convection speed for IU(x-tau*v).
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
  
  %Load Vector U(x-tauV) are obtained based on its position on the triangle
  %and a linear interpolation of 3 nodes of the triangle
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    flags = Mesh.ElemFlag; 
  else
    flags = zeros(nElements,1);
  end
  
     Mesh.ElemFlag = flags;
 
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
    
    %construct load corresponding to previous time step for IU(x-tau*v)
       
    [WU] = computeULagrangian_ConstantVelocity(U0,Vertices,Mesh,tau,Convectionhandle,Dhandle,StoreElement(i,:));
       
    
    L(idx(1)) = WU(1,1);
    L(idx(2)) = WU(2,1);
    L(idx(3)) = WU(3,1);
    
  end 