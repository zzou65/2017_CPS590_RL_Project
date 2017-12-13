function [Aloc] = STIMA_minSurfDer(Vertices,ElemInfo,w,varargin)
% STIMA_MINSURFDER Derivative matrix for the minimal surface problem when
% using a Newton scheme
%
%   ALOC = STIMA_MINSURFDER(VERTICES,ELEMINFO,W) computes the derivative
%   to element stiffness matrix with W as an estimate for solution U
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation
%
%   ELEMINFO is additional info about the elements
%
%   W is an estimate of the soulution U
%
%   Example:
%
%   Aloc = STIMA_minSurfDer([0 0 ; 1 0 ; 0 1],1,[1 2 3])

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland



  %Preallocate memory

  Aloc = zeros(3,3);

  %Compute the area of the element

  BK = [Vertices(2,:)-Vertices(1,:) ; ...
        Vertices(3,:)-Vertices(1,:)];
  Area = 0.5*abs(det(BK));

  inv_BK_t = transpose(inv(BK));

  % Compute element shape functions
  
  grad_N = grad_shap_LFE([0 0]);
  grad_N(1:2) = grad_N(1:2)*inv_BK_t;
  grad_N(3:4) = grad_N(3:4)*inv_BK_t;
  grad_N(5:6) = grad_N(5:6)*inv_BK_t;

  % Comoute entries of the stiffness matrix
  
  grad_U = w(1)*grad_N(1:2) + ...
           w(2)*grad_N(3:4) + ...
           w(3)*grad_N(5:6);
     
  z1 = sqrt(1 + sum(grad_U.^2));
  z2 = z1^3;

  Aloc(1,1) = (1/z1 * (grad_N(1:2) * grad_N(1:2)') - ...
               1/z2 * (grad_U * grad_N(1:2)') * (grad_U * grad_N(1:2)')) * Area;
  Aloc(1,2) = (1/z1 * (grad_N(1:2) * grad_N(3:4)') - ...
               1/z2 * (grad_U * grad_N(1:2)') * (grad_U * grad_N(3:4)')) * Area;
  Aloc(1,3) = (1/z1 * (grad_N(1:2) * grad_N(5:6)') - ...
               1/z2 * (grad_U * grad_N(1:2)') * (grad_U * grad_N(5:6)')) * Area;

  Aloc(2,2) = (1/z1 * (grad_N(3:4) * grad_N(3:4)') - ...
               1/z2 * (grad_U * grad_N(3:4)') * (grad_U * grad_N(3:4)')) * Area;
  Aloc(2,3) = (1/z1 * (grad_N(3:4) * grad_N(5:6)') - ...
               1/z2 * (grad_U * grad_N(3:4)') * (grad_U * grad_N(5:6)')) * Area;

  Aloc(3,3) = (1/z1 * (grad_N(5:6) * grad_N(5:6)') - ...
               1/z2 * (grad_U * grad_N(5:6)') * (grad_U * grad_N(5:6)')) * Area;

  % Fill in lower triangular part

  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);

return