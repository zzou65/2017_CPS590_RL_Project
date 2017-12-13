function Mloc = MASS_LFE2(Vertices,varargin)
% MASS_LFE2 Element mass matrix.
%
%   MLOC = MASS_LFE2(VERTICES) computes the element mass matrix using 
%   LFE2 finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_LFE2(Vertices);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  BK = [Vertices(2,:)-Vertices(1,:); (Vertices(3,:)-Vertices(1,:))];
  det_BK = abs(det(BK));

  % Compute local mass matrix
  
  Mloc = det_BK/24*[2 0 1 0 1 0;...
                    0 2 0 1 0 1;...
                    1 0 2 0 1 0;...
                    0 1 0 2 0 1;...
                    1 0 1 0 2 0;...
                    0 1 0 1 0 2];
   
return

