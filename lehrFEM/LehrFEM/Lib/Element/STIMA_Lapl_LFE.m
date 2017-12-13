function Aloc = STIMA_Lapl_LFE(Vertices,varargin)
% STIMA_LAPL_LFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_LAPL_LFE(VERTICES) computes the element stiffness matrix
%   for the Laplacian using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Aloc = STIMA_Lapl_LFE([0 0; 1 0; 0 1]);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(3,3);

  % Analytic computation of matrix entries using barycentric coordinates
  
  a = norm(Vertices(3,:)-Vertices(2,:)); 
  b = norm(Vertices(3,:)-Vertices(1,:));
  c = norm(Vertices(2,:)-Vertices(1,:));
  s = (a+b+c)/2;
  r = sqrt((s-a)*(s-b)*(s-c)/s);
  cot_1 = cot(2*atan(r/(s-a)));
  cot_2 = cot(2*atan(r/(s-b)));
  cot_3 = cot(2*atan(r/(s-c)));
  
  Aloc(1,1) = 1/2*(cot_3+cot_2);
  Aloc(1,2) = 1/2*(-cot_3);
  Aloc(1,3) = 1/2*(-cot_2);
  Aloc(2,2) = 1/2*(cot_3+cot_1);
  Aloc(2,3) = 1/2*(-cot_1);
  Aloc(3,3) = 1/2*(cot_2+cot_1);
  
  % Update lower triangular part
  
  Aloc(2,1) = Aloc(1,2);
  Aloc(3,1) = Aloc(1,3);
  Aloc(3,2) = Aloc(2,3);

return