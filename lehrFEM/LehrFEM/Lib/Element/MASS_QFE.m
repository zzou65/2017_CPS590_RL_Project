function Mloc = MASS_QFE(Vertices,varargin)
% MASS_QFE Element mass matrix.
%
%   MLOC = MASS_QFE(VERTICES) computes the element mass matrix using
%   quadratic Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_QFE(Vertices);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping
  
  BK = [Vertices(2,:)-Vertices(1,:); Vertices(3,:)-Vertices(1,:)];
  det_BK = abs(det(BK));
  
  % Fill in upper entries
  
  Mloc(1,1) = 1/60*det_BK;
  Mloc(1,2) = -1/360*det_BK;
  Mloc(1,3) = -1/360*det_BK ;
  Mloc(1,4) = 0;
  Mloc(1,5) = -1/90*det_BK;
  Mloc(1,6) = 0;
  
  Mloc(2,2) = 1/60*det_BK;
  Mloc(2,3) = -1/360*det_BK;
  Mloc(2,4) = 0;
  Mloc(2,5) = 0;
  Mloc(2,6) = -1/90*det_BK;
  
  Mloc(3,3) = 1/60*det_BK;
  Mloc(3,4) = -1/90*det_BK;
  Mloc(3,5) = 0;
  Mloc(3,6) = 0;
  
  Mloc(4,4) = 4/45*det_BK;
  Mloc(4,5) = 2/45*det_BK;
  Mloc(4,6) = 2/45*det_BK;
 
  Mloc(5,5) = 4/45*det_BK;
  Mloc(5,6) = 2/45*det_BK;
 
  Mloc(6,6) = 4/45*det_BK;
  
  % Fill in lower triangular part
  
  tri = triu(Mloc);
  Mloc = tri+tril(tri',-1);
 
return