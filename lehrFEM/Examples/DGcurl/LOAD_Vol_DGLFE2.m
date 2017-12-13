function Lloc = LOAD_Vol_DGLFE2(Vertices,ElemFlag,QuadRule,FHandle1,FHandle2,varargin)
% LOAD_VOL_DGLFE2 Element load vector for volume load data.
%
%   LLOC = LOAD_VOL_DGLFE2(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,S,FHANDLE)
%   computes the entries of the element load vector for the boundary load
%   data.
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


%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  
  nPts = size(QuadRule.x,1);

  % Preallocate memory
  
  Lloc = zeros(6,1);
 
  % Compute values of shape functions
  
  N = shap_DGLFE(QuadRule.x);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; ...
        Vertices(3,:)-bK];
  det_BK = abs(det(BK));
    
  x = QuadRule.x*BK+ones(nPts,1)*bK;
  
  
  % Compute function values
  
  FVal1 = FHandle1(x,ElemFlag,varargin{:});
  FVal2 = FHandle2(x,ElemFlag,varargin{:});
  FVal = [FVal1 FVal2];
  
  % Compute entries of element load vector
   
  Lloc(1) = sum(QuadRule.w.*FVal(:,1).*N(:,1))*det_BK;
  Lloc(2) = sum(QuadRule.w.*FVal(:,2).*N(:,1))*det_BK;  
  Lloc(3) = sum(QuadRule.w.*FVal(:,1).*N(:,2))*det_BK;
  Lloc(4) = sum(QuadRule.w.*FVal(:,2).*N(:,2))*det_BK;
  Lloc(5) = sum(QuadRule.w.*FVal(:,1).*N(:,3))*det_BK;
  Lloc(6) = sum(QuadRule.w.*FVal(:,2).*N(:,3))*det_BK;
  
return
