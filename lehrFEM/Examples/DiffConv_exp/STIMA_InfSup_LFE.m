function Aloc = STIMA_InfSup_LFE(Vertices, flag, QuadRule, VHandle, varargin)
% STIMA_SUPG_LFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_SUPG_LFE(VERTICES) computes the element stiffness matrix
%   for SUPG-modification for convection using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%  
%   a: diffusivity
%   d1 d2: apriori chosen constants for SUPG-modification
%
%   Flag useless, needed for interface to assemMat_LFE
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.e: 
%
%   VHANDLE is function handle for velocity field   
%
%   Example:
%
%   Aloc = STIMA_Lapl_LFE([0 0; 1 0; 0 1]);

%   Copyright 2005-2007 Patrick Meury, Holger Heumann, Alan Mitchell
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
% This stiffness matrix provides the required modification to the load matrix on the RHS for
% the SUPG method


  % Preallocate memory
  
  Aloc = zeros(3,3);

  % Analytic computation of matrix entries using barycentric coordinates
  
  l1x = Vertices(2,2)-Vertices(3,2); 
  l1y = Vertices(3,1)-Vertices(2,1);
  l2x = Vertices(3,2)-Vertices(1,2);
  l2y = Vertices(1,1)-Vertices(3,1); 
  l3x = Vertices(1,2)-Vertices(2,2);
  l3y = Vertices(2,1)-Vertices(1,1);
   
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ];          % transpose of transformation matrix
  det_BK = abs(det(BK));               % twice the area of the triagle
    
  nPoints = size(QuadRule.w,1);
  
  % Quadrature points in actual element stored as rows of a matrix
  x = QuadRule.x*BK + ones(nPoints,1)*P1;

  % Evaluate coefficient function at quadrature nodes
  c = VHandle(x,varargin{:});
  FHandle = [c(:,1).*c(:,1) c(:,1).*c(:,2) c(:,2).*c(:,1) c(:,2).*c(:,2)];
   
  % Apply quadrature rule and fix constant part
  
  w = QuadRule.w;
  e = sum((FHandle.*[w w w w]), 1);
  
  te(1,1) = e(1);
  te(1,2) = e(2);
  te(2,1) = e(3);
  te(2,2) = e(4);
  te = te./det_BK;
 
   % Compute Aloc values
  
  Aloc(1,1) = (te*[l1x l1y]')'*[l1x l1y]';            
  Aloc(1,2) = (te*[l1x l1y]')'*[l2x l2y]';            
  Aloc(1,3) = (te*[l1x l1y]')'*[l3x l3y]';            
  Aloc(2,2) = (te*[l2x l2y]')'*[l2x l2y]';            
  Aloc(2,3) = (te*[l2x l2y]')'*[l3x l3y]';            
  Aloc(3,3) = (te*[l3x l3y]')'*[l3x l3y]';            
  Aloc(2,1) = (te*[l2x l2y]')'*[l1x l1y]';            
  Aloc(3,1) = (te*[l3x l3y]')'*[l1x l1y]';            
  Aloc(3,2) = (te*[l3x l3y]')'*[l2x l2y]';
  
return