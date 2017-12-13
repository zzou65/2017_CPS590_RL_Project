function Aloc = STIMA_Curl_LFE2(Vertices,varargin)
% STIMA_CURL_LFE2 element stiffness matrix.
%
%   ALOC = STIMA_CURL_LFE2(VERTICES) computes the element stiffness matrix 
%   using nodal finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   Example:
%
%   Aloc = STIMA_Curl_LFE2(Vertices);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute the area of the element
  
  BK = [Vertices(2,:)-Vertices(1,:);Vertices(3,:)-Vertices(1,:)];
  det_BK = abs(det(BK));
  
  % Compute local mass matrix
  
  K = [ Vertices(3,:) - Vertices(2,:) ...
        Vertices(1,:) - Vertices(3,:) ...
        Vertices(2,:) - Vertices(1,:) ];
   
  Aloc = 1/(2*det_BK)*(K')*K;
   
return
