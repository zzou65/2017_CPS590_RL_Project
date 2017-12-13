function [M,S,L,storeElementMatrix] = SemiLagrangian_LFE(Mesh,U0,tau,Convectionhandle,Dhandle,QuadRule)
% SemiLagrangian_LFE Asseble Mass Matrix, Laplacian Matrix and Load Vector 
% Corresponding to previous time step and lagrangian tracking position.
%
%   L = SemiLagrangian_LFE(Mesh,U0,tau,Convectionhandle,Dhandle,QuadRule) assembles  
%   Mass Matrix, laplacian matrix and Load Vector based on Convectionhandle
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
  
  nPts = size(QuadRule.w,1);
  storeElementMatrix = zeros(nElements,nPts);
  
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

  N = shap_LFE(QuadRule.x); 
  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nElements
    
    % Extract vertices of current element
    
    idx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(idx,:);
    
    %construct load corresponding to previous time step
    
    %construct LFE using U and linear functions
       
    [WU,storeElement] = computeLoadIntegral(U0,QuadRule,Vertices,Mesh,tau,Convectionhandle,Dhandle);
    storeElementMatrix(i,:) = storeElement;
    
    bK = Mesh.Coordinates(idx(1),:);
    BK = [Mesh.Coordinates(idx(2),:)-bK; Mesh.Coordinates(idx(3),:)-bK];
    det_BK = abs(det(BK));
      
    
    L(idx(1)) = L(idx(1)) + sum(WU.*N(:,1))*det_BK;
    L(idx(2)) = L(idx(2)) + sum(WU.*N(:,2))*det_BK;
    L(idx(3)) = L(idx(3)) + sum(WU.*N(:,3))*det_BK;
    
    
    % load is completed
      
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
    