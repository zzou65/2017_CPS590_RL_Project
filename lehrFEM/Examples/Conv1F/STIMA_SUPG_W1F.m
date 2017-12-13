function Mloc = STIMA_SUPG_W1F(Vertices,ElemInfo,V_HANDLE,QuadRule,varargin)
% STIMA_SUPG_W1F 
%
%   MLOC = STIMA_SUPG_W1F(VERTICES ...) computes ... matrix using 
%   Whitney 1-forms finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
% 
%   ElemInfo (not used)
%
%   V_HANDLE 
%   Example:
%

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ]; % transpose of transformation matrix
  det_BK = abs(det(BK));     % twice the area of the triagle
  nPoints=size(QuadRule.x,1);
  
  x = QuadRule.x*BK + ones(nPoints,1)*P1;

  % Evaluate coefficient function at quadrature nodes
  Fval = V_HANDLE(x,ElemInfo,varargin{:});
  
    
  % Compute local mass matrix
  weights = QuadRule.w*4/det_BK;% * det_BK/det_Bk^2*4;
  Mloc=sum(sum(Fval.*Fval,2).*weights)*ones(3,3);
  return