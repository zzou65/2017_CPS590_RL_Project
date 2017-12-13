function M = MASS_Weight_LFE(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% MASS_WEIGHT_LFE Element mass matrix.
%
%   M = MASS_WEIGHT_LFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the 
%   element mass matrix using bi-linear Lagrangian finite elements.
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used to
%   do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   FHANDLE denotes the function handle to the weight functions.
%
%   M = MASS_WEIGHT_LFE(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also handles
%   the additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   M = MASS_Weight_LFE([0 0; 1 0; 0 1],0,P7O6(),FHandle);
%
%   See also MASS_BFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);

  % Preallocate memory
  
  M = zeros(3,3);
  
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  N = shap_LFE(QuadRule.x);

  bK = P1;
  BK = [ P2 - P1 ; P3 - P1 ];
  det_BK = abs(det(BK));
  
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  Fval = FHandle(x,ElemInfo,varargin{:});
  
  % Compute the matrix
  
  for i = 1:3
    for j = i:3
      M(i,j) = sum(QuadRule.w.*Fval.*N(:,i).*N(:,j))*det_BK; 
    end
  end
  
  % Fill in lower triangular part
  
  tri = triu(M);
  M = tri+tril(tri.',-1);
  
  
return