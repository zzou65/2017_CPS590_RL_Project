function Mloc = MASS_LFE(Vertices,varargin)
% MASS_LFE Element mass matrix.
%
%   MLOC = MASS_LFE(VERTICES) computes the element mass matrix using linear
%   Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_LFE(Vertices);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; ...
        Vertices(3,:)-bK];
  det_BK = abs(det(BK));

  % Compute local mass matrix
  
  Mloc = det_BK/24*[2 1 1; 1 2 1; 1 1 2];
   
return