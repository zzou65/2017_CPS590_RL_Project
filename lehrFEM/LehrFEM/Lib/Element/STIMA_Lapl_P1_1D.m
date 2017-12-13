function Aloc = STIMA_Lapl_P1_1D(Vertices,varargin)
% STIMA_LAPL_LFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_LAPL_P1_1D(VERTICES) computes the element stiffness matrix
%   for the Laplacian using linear Lagrangian finite elements.
%
%   VERTICES is 2-by-1 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Aloc = STIMA_Lapl_P1_1D([0; 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  h = abs(Vertices(2)-Vertices(1));
  
  Aloc = 1/h*ones(2,2);
  Aloc(1,2) = -Aloc(1,2);
  Aloc(2,1) = -Aloc(2,1);
  
return