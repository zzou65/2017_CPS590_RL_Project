function Mloc = MASS_Vol_DGLFE(Vertices,varargin)
% MASS_VOL_DGLFE Element mass matrix.
%
%   MLOC = MASS_VOL_DGLFE(VERTICES) computes the element mass matrix using
%   discontinuous linear finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_Vol_DGLFE(Vertices);

%   Copyright 2006-2006 Patrick Meury
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