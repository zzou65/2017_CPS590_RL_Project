function Aloc = MASS_Lump_LFE(Vertices,varargin)
% MASS_LUMP_LFE element lumbda mass matrix.
%
%   ALOC = MASS_LUMP_LFE(VERTICES) computes the element mass matrix 
%   using W1F finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   Example:
%
%   Aloc = MASS_Lump_LFE(Vertices);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute the area of the element
  
  BK = [Vertices(2,:)-Vertices(1,:);Vertices(3,:)-Vertices(1,:)];
  det_BK = abs(det(BK));
  
  % Compute local mass matrix
  
  Aloc = 1/6*det_BK*eye(3);
   
return