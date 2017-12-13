function Lloc = LOAD_Vol_DGCR(Vertices,ElemFlag,QuadRule,FHandle,varargin)
% LOAD_VOL_DG Element load vector for volume load data.
%
%   LLOC = LOAD_VOL_DGCR(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,S,FHANDLE) computes
%   the entries of the element load vector for the boundary load data.
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
%   FHANDLE is a function pointer to the load data.
%
%   LLOC = LOAD_VOL_DGCR(VERTICES,ELEMFLAG,QUADRULE,FHANDLE,FPARAM) also
%   handles the variable length argumet list FPARAM to the function pointer
%   DHANDLE.
%
%   Example:
%
%   F = @(x,varargin)-4*ones(size(x,1),1);
%   Lloc = LOAD_Vol_DGCR(Vertices,ElemFlag,P3O3(),F);
%
%   See also shap_DGCR.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  
  nPts = size(QuadRule.x,1);

  % Preallocate memory
  
  Lloc = zeros(3,1);
 
  % Compute values of shape functions
  
  N = shap_DGCR(QuadRule.x);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; ...
        Vertices(3,:)-bK];
  det_BK = abs(det(BK));
    
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  % Compute function values
  
  FVal = FHandle(x,ElemFlag,varargin{:});
  
  % Compute entries of element load vector
   
  Lloc(1) = sum(QuadRule.w.*FVal.*N(:,1))*det_BK;
  Lloc(2) = sum(QuadRule.w.*FVal.*N(:,2))*det_BK;
  Lloc(3) = sum(QuadRule.w.*FVal.*N(:,3))*det_BK;
  
return
