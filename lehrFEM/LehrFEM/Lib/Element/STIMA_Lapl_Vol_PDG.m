function Aloc = STIMA_Lapl_Vol_PDG(Vertices,ElemFlag,QuadRule,grad_Shap,varargin)
% STIMA_LAPL_VOL_PDG Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_LAPL_VOL_PDG(VERTICES,ELEMFLAG,QUADRULE,GRAD_SHAP) computes
%   the element stiffness matrix for the Laplacian using the shape functions
%   given by the function handle GRAD_SHAP.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   The integer ELEMFLAG denotes the element flag of the current element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   GRAD_SHAP is a function handle to the gradients of the reference
%   element shape functions.
%
%   Example:
%
%   Aloc = STIMA_Lapl_Vol_DGCR(Vertices,ElemFlag,P3O3(),@grad_shap_DGCR);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Compute element mapping

  bK = Vertices(1,:);
  BK = [Vertices(2,:) - bK; ...
        Vertices(3,:) - bK];
  det_BK = abs(det(BK));
  inv_BK = inv(BK);
  
  TK = det_BK*transpose(inv_BK)*inv_BK;
  
  % Evaluate shape functions
  
  grad_N = grad_Shap(QuadRule.x);
  nDofs = size(grad_N,2)/2;
  
  % Preallocate memory
  
  Aloc = zeros(nDofs,nDofs);
  
  % Compute local mass matrix
  
  for j1 = 1:nDofs
    loc_1 = 2*(j1-1) + [1 2];
    for j2 = j1:nDofs
      loc_2 = 2*(j2-1) + [1 2];
      Aloc(j1,j2) = sum(QuadRule.w.*sum(grad_N(:,loc_1).*(grad_N(:,loc_2)*TK),2));
    end
  end

  % Fill in lower trinagular part
    
  tri = triu(Aloc);
  Aloc = tril(transpose(tri),-1)+tri;
  
return