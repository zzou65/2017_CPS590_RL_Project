function Aloc = STIMA_minSurf(Vertices,ElemInfo,w,varargin)
% STIMA_MINSURF Elements stiffness matrix for the minimal surface problem
%
%   ALOC = STIMA_MINSURF(VERTICES,ELEMINFO,W) computes the element stiffness matrix
%   with W as an estimate for solution U
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation
%
%   ELEMINFO is additional info about the elements
%
%   W is an estimate of the solution U
%
%   Example:
%
%   Aloc = STIMA_minSurf([0 0 ; 1 0 ; 0 1],1,[1 2 3])

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Preallocate memory
  
  Aloc = zeros(3,3);
  
  % Compute the area of the element

  BK = [Vertices(2,:)-Vertices(1,:); ...
        Vertices(3,:)-Vertices(1,:)];
  Area = 0.5*abs(det(BK));
    
  inv_BK_t = transpose(inv(BK));
  
  % Compute element shape functions

  grad_N = grad_shap_LFE([0 0]);
  grad_N(1:2) = grad_N(1:2)*inv_BK_t;
  grad_N(3:4) = grad_N(3:4)*inv_BK_t;
  grad_N(5:6) = grad_N(5:6)*inv_BK_t;
  
  % Compute entries of the stfiffness matrix
  
  z = sqrt(1 + sum((w(1)*grad_N(1:2) + ...
                    w(2)*grad_N(3:4) + ...
                    w(3)*grad_N(5:6)).^2));
  
  Aloc(1,1) = sum(grad_N(1:2).*grad_N(1:2))*Area/z;
  Aloc(1,2) = sum(grad_N(1:2).*grad_N(3:4))*Area/z;
  Aloc(1,3) = sum(grad_N(1:2).*grad_N(5:6))*Area/z;
  
  Aloc(2,2) = sum(grad_N(3:4).*grad_N(3:4))*Area/z;
  Aloc(2,3) = sum(grad_N(3:4).*grad_N(5:6))*Area/z;
  
  Aloc(3,3) = sum(grad_N(5:6).*grad_N(5:6))*Area/z;
  
  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
  
return