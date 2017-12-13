function M = MASS_P0(Vertices,ElemInfo,QuadRule,FHandle,varargin)
% MASS_P0 Element mass matrix.
%
%   M = MASS_P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE) computes the element
%   mass matrix using constant Lagrangian finite elements.
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   M = MASS_P0(VERTICES,ELEMINFO,QUADRULE,FHANDLE,FPARAM) also handles the
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   M = MASS_P0([0 0; 1 0; 0 1],0,P7O6,FHandle);
%
%   See also MASS_P0.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK ; Vertices(3,:)-bK];
  det_BK = abs(det(BK));
  
  % Compute the matrix
  
  M = det_BK/2; 
  
return