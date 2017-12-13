function Mloc = MASSlump_P1_1D(Vertices,varargin)
% MASS_P1_1D Element mass matrix.
%
%   MLOC = MASS_P1_1D(VERTICES) computes the lumped element mass matrix using
%   linear finite elements.
%
%   VERTICES is 2-by-1 matrix specifying the vertices of the current
%   element.
%
%   Example:
%
%   Mloc = MASS_P1_1D([0; 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  h = abs(Vertices(2)-Vertices(1));
  
  Mloc = zeros(2,2);
  Mloc(1,1) = h/2;
  Mloc(2,2) = h/2;

return
