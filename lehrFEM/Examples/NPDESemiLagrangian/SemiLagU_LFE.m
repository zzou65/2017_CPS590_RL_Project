function [M,S,L,storeElementMatrix] = SemiLagU_LFE(Mesh,U0,tau,Convectionhandle,Dhandle)
% SemiLagU_LFE Asseble Mass Matrix, Laplacian Matrix and L Vector 
% Corresponding to previous time step and lagrangian tracking position for IU(x-tau*v).
%
%   [M,S,L,storeElementMatrix] = SemiLagU_LFE(Mesh,U0,tau,Convectionhandle,Dhandle,QuadRule) assembles  
%   Mass Matrix, laplacian matrix and Load Vector based on Convectionhandle
%   and storeElementMatrix to be used for time independent velocity field
%   in later time steps
%
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
  Mhandle = @MASS_LFE;
  Shandle = @STIMA_Lapl_LFE;
  
  %Load Vector U(x-tauV)*V(x) are integrated based on P3O2() Gaussian
  %Quadrature
  
  storeElementMatrix = zeros(nElements,3);
  
  % Preallocate memory
  
  MI = zeros(9*nElements,1);
  MJ = zeros(9*nElements,1);
  MA = zeros(9*nElements,1);
  
  SI = zeros(9*nElements,1);
  SJ = zeros(9*nElements,1);
  SA = zeros(9*nElements,1);
  L = zeros(nCoordinates,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    flags = Mesh.ElemFlag; 
  else
    flags = zeros(nElements,1);
  end
  
     Mesh.ElemFlag = flags;

  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
    
    %construct load corresponding to previous time step for IU(x-tau*v)
       
    [WU,storeElement] = computeULagrangian(U0,Vertices,Mesh,tau,Convectionhandle,Dhandle);
    storeElementMatrix(i,:) = storeElement;
        
    L(idx(1)) = WU(1,1);
    L(idx(2)) = WU(2,1);
    L(idx(3)) = WU(3,1);
      
    % Compute element contributions
    
    Mloc = Mhandle(Vertices,flags(i));
    Sloc = tau*Shandle(Vertices,flags(i));
   
    % Add contributions to stiffness matrix
    
    MI(loc) = set_Rows(idx,3);
    MJ(loc) = set_Cols(idx,3);
    MA(loc) = Mloc(:);
    
    SI(loc) = set_Rows(idx,3);
    SJ(loc) = set_Cols(idx,3);
    SA(loc) = Sloc(:);
    
    loc = loc+9;
    
  end
  
  % Assign output arguments
    M = sparse(MI,MJ,MA);    
    
    S = sparse(SI,SJ,SA); 