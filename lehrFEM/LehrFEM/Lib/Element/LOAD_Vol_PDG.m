function Lloc = LOAD_Vol_PDG(Vertices,ElemFlag,QuadRule,Shap,FHandle,varargin)
% LOAD_VOL_PDG Element load vector for volume load data.
%
%   LLOC = LOAD_VOL_PDG(VERTICES,ELEMFLAG,QUADRULE,SHAP,FHANDLE) computes
%   the volume contributions of the element load vector using the shape 
%   functions given by the function handle SHAP.
%
%   VERTICES is 3-by-2 matrix whose rows contain the vertices of the
%   current element.
%
%   The integer ELEMFLAG denotes the element flag of the current element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHAP is a function handle to the reference element shape functions.
%
%   FHANDLE is a function handle to the load data.
%
%   LLOC = LOAD_VOL_PDG(VERTICES,ELEMFLAG,QUADRULE,SHAP,FHANDLE,FPARAM)
%   also handles the variable length argumet list FPARAM to the function
%   pointer FHANDLE.
%
%   Example:
%
%   F = @(x,varargin)-4*ones(size(x,1),1);
%   Lloc = LOAD_Vol_PDG(Vertices,ElemFlag,P3O3(),@shap_DGCR,F);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  
  nPts = size(QuadRule.x,1);  % Number of quadrature points
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:) - bK; ...
        Vertices(3,:) - bK];
  det_BK = abs(det(BK));
    
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Evaluate shape functions
  
  N = Shap(QuadRule.x);
  nDofs = size(N,2);
  
  % Preallocate memory
  
  Lloc = zeros(nDofs,1);
  
  % Compute function values
  
  FVal = FHandle(x,ElemFlag,varargin{:});
  
  % Compute entries of element load vector
  
  for j = 1:nDofs
    Lloc(j) = sum(QuadRule.w.*FVal.*N(:,j))*det_BK;
  end
  
return
