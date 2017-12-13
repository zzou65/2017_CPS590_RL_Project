function Mloc = MASS_P0_1D(Vertices,varargin)
% MASS_P0_1D Element mass matrix.
%
%   MLOC = MASS_P0_1D(VERTICES) computes the element mass matrix using
%   linear finite elements.
%
%   VERTICES is 2-by-1 matrix specifying the vertices of the current
%   element.
%
%   Example:
%
%   Mloc = MASS_P0_1D([0; 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  h = abs(Vertices(2)-Vertices(1));
  Mloc = h;
  
return