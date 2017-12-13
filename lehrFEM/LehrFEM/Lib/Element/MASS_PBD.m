function Mloc = MASS_PBD(Vertices,ElemFlag,delta,QuadRule,varargin)
% MASS_PBD Element mass matrix.
%
%   MLOC = MASS_PBD(VERTICES,ELEMFLAG,DELTA,QUADRULE) computes the element
%   mass matrix using quadratic Lagrangian finite elements together with
%   parabolic boundary corrections.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%   
%   ELEMFLAG is an integer parameter which is used to specify additional
%   element information on each element.
%
%   DELTA specifies the value of the parabolic boundary correction term on
%   the first edge.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   Example:
%
%   Mloc = MASS_PBD([0 0; 1 0; 0 1],0,0.25,P7O6());
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  Rot = [0 -1; 1 0];
  nPts = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Mloc = zeros(6,6);
  
  % Precompute shape fuction values
  
  N = shap_QFE(QuadRule.x);
    
  % Compute standard element mapping
 
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
   
  if(delta > eps)
  
    % Compute element mass matrix for curved elements
    
    normal = Vertices(3,:)-Vertices(2,:);
    normal = normal*Rot/norm(normal);
    for i = 1:nPts      
      Jac = transpose(BK) + 4*delta*transpose(normal) * ...
            [-QuadRule.x(i,2) 1-QuadRule.x(i,1)-2*QuadRule.x(i,2)];  
      det_Jac = abs(det(Jac));
      for j1 = 1:6
        for j2 = j1:6  
          Mloc(j1,j2) = Mloc(j1,j2) + QuadRule.w(i)*N(i,j1)*N(i,j2)*det_Jac;
        end
      end
    end
    
  else
    
    % Compute element mass matrix for straight elements
      
    det_BK = abs(det(BK));     
    for j1 = 1:6
      for j2 = j1:6
        Mloc(j1,j2) = sum(QuadRule.w.*N(:,j1).*N(:,j2))*det_BK;  
      end
    end    
        
  end
  
  % Fill in lower triangular part
  
  tri = triu(Mloc);
  Mloc = tri+tril(tri',-1);
  
return