function Aloc = MASS_Lump_Inv_LFE(Vertices,varargin)
% MASS_LUMP_LFE element lump mass matrix.
%
%   ALOC = MASS_LUMP_Inv_LFE(VERTICES) computes the inverse of element mass matrix 
%   using nodal finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute the area of the element
  
  BK = [Vertices(2,:)-Vertices(1,:);Vertices(3,:)-Vertices(1,:)];
  det_BK = abs(det(BK));
  
  % Compute local mass matrix
  
  Aloc = 6/(det_BK)*eye(3);
   
return