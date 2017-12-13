function Aloc = STIMA_Lapl_QFE(Vertices,ElemInfo,varargin)
% STIMA_LAPL_QFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_LAPL_QFE(VERTICES,ELEMINFO) computes the element stiffness
%   matrix for the Laplacian using quadratic Lagrangian finite elements.
%
%   VERTICES is 3-by-3 matrix specifying the vertices of the current 
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   Example:
%
%   Aloc = STIMA_Lapl_LFE(Vertices,ElemInfo);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize local stiffness matrix
   
  Aloc = zeros(6,6);
  
  % Analytic computation of matrix entries using barycentric coordinates
  
  a = norm(Vertices(3,:)-Vertices(2,:)); 
  b = norm(Vertices(3,:)-Vertices(1,:));
  c = norm(Vertices(2,:)-Vertices(1,:));
  s = (a+b+c)/2;
  r = sqrt((s-a)*(s-b)*(s-c)/s);
  cot_1 = cot(2*atan(r/(s-a)));
  cot_2 = cot(2*atan(r/(s-b)));
  cot_3 = cot(2*atan(r/(s-c)));
  
  Aloc(1,1) = (1/2)*(cot_3+cot_2);
  Aloc(1,2) = (-1/6)*(-cot_3);
  Aloc(1,3) = (-1/6)*(-cot_2);
  Aloc(1,4) = (2/3)*(-cot_3);
  Aloc(1,5) = 0;
  Aloc(1,6) = (2/3)*(-cot_2);
    
  Aloc(2,2) = 1/2*(cot_3+cot_1);
  Aloc(2,3) = (-1/6)*(-cot_1);
  Aloc(2,4) = (2/3)*(-cot_3);
  Aloc(2,5) = (2/3)*(-cot_1);
  Aloc(2,6) = 0;
  
  Aloc(3,3) = 1/2*(cot_2+cot_1);
  Aloc(3,4) = 0;
  Aloc(3,5) = (2/3)*(-cot_1);
  Aloc(3,6) = (2/3)*(-cot_2);
  
  Aloc(4,4) = (4/3)*(cot_3+cot_2+cot_1);
  Aloc(4,5) = (4/3)*(-cot_2);
  Aloc(4,6) = (4/3)*(-cot_1);
  
  Aloc(5,5) = (4/3)*(cot_3+cot_2+cot_1);
  Aloc(5,6) = (4/3)*(-cot_3);
  
  Aloc(6,6) = (4/3)*(cot_3+cot_2+cot_1);
  
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return