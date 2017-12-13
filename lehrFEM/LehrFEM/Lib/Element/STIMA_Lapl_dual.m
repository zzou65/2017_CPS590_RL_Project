function Aloc = STIMA_Lapl_dual(Vertices,ElemInfo,QuadRule,varargin)
% STIMA_Lapl_dual Element stiffness matrix for the laplace problem in dual formulation.
%
%   ALOC = STIMA_Lapl_dual(VERTICES,ELEMINFO,QUADRULE) computes the
%   element stiffness matrix for the laplace problem in dual formulation. 
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   Example:
%
%   Aloc = STIMA_Lapl_dual([0 0; 1 0; 0 1],0,P7O6());
%

%   Copyright 2005-2006 Patrick Meury & Kah Ling Sia & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  MU_Handle = @(x,varargin)ones(size(x,1),1);
  Aloc = zeros(4,4);
  
  % Compute element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK_t = transpose(inv(BK));
  det_BK = abs(det(BK));

  % Compute element stiffness matrix
  
  % mass term of flux discretisation
  Aloc(1:3,1:3)=MASS_W1F(Vertices,0,MU_Handle,QuadRule,varargin);
  
  % upper off diagonal block
   
  Aloc(4,1)=1;
  Aloc(4,2)=1;
  Aloc(4,3)=1;
  
  % lower off diagonal block
  
  Aloc(1,4)=1;
  Aloc(2,4)=1;
  Aloc(3,4)=1;
return