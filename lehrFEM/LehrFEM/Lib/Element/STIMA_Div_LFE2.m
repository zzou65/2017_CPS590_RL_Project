function Aloc = STIMA_Div_LFE2(Vertices,varargin)
% STIMA_DIV_LFE2 element stiffness matrix.
%
%   ALOC = STIMA_DIV_LFE2(VERTICES,ELEMINFO,F_HANDLE,QUADRULE) computes the
%   element stiffness matrix using nodal finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Aloc = STIMA_Div_LFE2(Vertices);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute the area of the element
  
  BK = [Vertices(2,:)-Vertices(1,:);Vertices(3,:)-Vertices(1,:)];
  det_BK = abs(det(BK));
  
  % Compute local mass matrix
  
  L = [ Vertices(2,2) - Vertices(3,2) ...
        Vertices(3,1) - Vertices(2,1) ...
        Vertices(3,2) - Vertices(1,2) ...
        Vertices(1,1) - Vertices(3,1) ...
        Vertices(1,2) - Vertices(2,2) ...
        Vertices(2,1) - Vertices(1,1) ];
   
  Aloc = 1/(2*det_BK)*L'*L;
   
return
