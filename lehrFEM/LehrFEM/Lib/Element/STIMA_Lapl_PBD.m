function Aloc = STIMA_Lapl_PBD(Vertices,ElemFlag,delta,QuadRule,varargin)
% STIMA_LAPL_PBD Element mass matrix.
%
%   MLOC = STIMA_LAPL_PBD(VERTICES,ELEMFLAG,DELTA,QUADRULE) computes the
%   element mass matrix using linear Lagrangian finite elements together
%   with parabolic boundary corrections.
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
%   Mloc = STIMA_Lapl_PBD([0 0; 1 0; 0 1],0,0.25,P7O6());
%
%   See also grad_shap_QFE.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  Rot = [0 -1; 1 0];
  nPts = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Aloc = zeros(6,6);
  
  % Precompute shape function values
  
  grad_N = grad_shap_QFE(QuadRule.x);
    
  % Compute standard element mapping
 
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
   
  if(delta > eps)
  
    % Compute element stiffness matrix for curved elements
    
    normal = Vertices(3,:)-Vertices(2,:);
    normal = normal*Rot/norm(normal);
    for i = 1:nPts      
      Jac = transpose(BK) + 4*delta*transpose(normal) * ...
            [-QuadRule.x(i,2) 1-QuadRule.x(i,1)-2*QuadRule.x(i,2)];  
      inv_Jac = inv(Jac);
      TK = abs(det(Jac))*inv_Jac*transpose(inv_Jac);
      
      Aloc(1,1) = Aloc(1,1) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,1:2));
      Aloc(1,2) = Aloc(1,2) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,3:4));
      Aloc(1,3) = Aloc(1,3) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,5:6));
      Aloc(1,4) = Aloc(1,4) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,7:8));
      Aloc(1,5) = Aloc(1,5) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,9:10));
      Aloc(1,6) = Aloc(1,6) + ...
                  QuadRule.w(i)*grad_N(i,1:2)*TK*transpose(grad_N(i,11:12));
      
       Aloc(2,2) = Aloc(2,2) + ...
                   QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,3:4));       
       Aloc(2,3) = Aloc(2,3) + ...
                   QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,5:6));
       Aloc(2,4) = Aloc(2,4) + ...
                   QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,7:8));
       Aloc(2,5) = Aloc(2,5) + ...
                   QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,9:10));
       Aloc(2,6) = Aloc(2,6) + ...
                   QuadRule.w(i)*grad_N(i,3:4)*TK*transpose(grad_N(i,11:12));
              
       Aloc(3,3) = Aloc(3,3) + ...
                   QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,5:6));
       Aloc(3,4) = Aloc(3,4) + ...
                   QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,7:8));
       Aloc(3,5) = Aloc(3,5) + ...
                   QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,9:10));
       Aloc(3,6) = Aloc(3,6) + ...
                   QuadRule.w(i)*grad_N(i,5:6)*TK*transpose(grad_N(i,11:12));
              
       Aloc(4,4) = Aloc(4,4) + ...
                   QuadRule.w(i)*grad_N(i,7:8)*TK*transpose(grad_N(i,7:8));
       Aloc(4,5) = Aloc(4,5) + ...
                   QuadRule.w(i)*grad_N(i,7:8)*TK*transpose(grad_N(i,9:10));
       Aloc(4,6) = Aloc(4,6) + ...
                   QuadRule.w(i)*grad_N(i,7:8)*TK*transpose(grad_N(i,11:12));
              
       Aloc(5,5) = Aloc(5,5) + ...
                   QuadRule.w(i)*grad_N(i,9:10)*TK*transpose(grad_N(i,9:10));
       Aloc(5,6) = Aloc(5,6) + ... 
                   QuadRule.w(i)*grad_N(i,9:10)*TK*transpose(grad_N(i,11:12));
              
       Aloc(6,6) = Aloc(6,6) + ...
                   QuadRule.w(i)*grad_N(i,11:12)*TK*transpose(grad_N(i,11:12));
               
    end
    
  else
        
    % Compute element stiffness matrix for straight elements
      
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    TK = det_BK*transpose(inv_BK)*inv_BK;
    
    Aloc(1,1) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,1:2)*TK),2));
    Aloc(1,2) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,3:4)*TK),2));
    Aloc(1,3) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,5:6)*TK),2));
    Aloc(1,4) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,7:8)*TK),2));
    Aloc(1,5) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,9:10)*TK),2));
    Aloc(1,6) = sum(QuadRule.w.*sum((grad_N(:,1:2)).*(grad_N(:,11:12)*TK),2));
  
    Aloc(2,2) = sum(QuadRule.w.*sum((grad_N(:,3:4)).*(grad_N(:,3:4)*TK),2));
    Aloc(2,3) = sum(QuadRule.w.*sum((grad_N(:,3:4)).*(grad_N(:,5:6)*TK),2));
    Aloc(2,4) = sum(QuadRule.w.*sum((grad_N(:,3:4)).*(grad_N(:,7:8)*TK),2));
    Aloc(2,5) = sum(QuadRule.w.*sum((grad_N(:,3:4)).*(grad_N(:,9:10)*TK),2));
    Aloc(2,6) = sum(QuadRule.w.*sum((grad_N(:,3:4)).*(grad_N(:,11:12)*TK),2));
  
    Aloc(3,3) = sum(QuadRule.w.*sum((grad_N(:,5:6)).*(grad_N(:,5:6)*TK),2));
    Aloc(3,4) = sum(QuadRule.w.*sum((grad_N(:,5:6)).*(grad_N(:,7:8)*TK),2));
    Aloc(3,5) = sum(QuadRule.w.*sum((grad_N(:,5:6)).*(grad_N(:,9:10)*TK),2));
    Aloc(3,6) = sum(QuadRule.w.*sum((grad_N(:,5:6)).*(grad_N(:,11:12)*TK),2));
  
    Aloc(4,4) = sum(QuadRule.w.*sum((grad_N(:,7:8)).*(grad_N(:,7:8)*TK),2));
    Aloc(4,5) = sum(QuadRule.w.*sum((grad_N(:,7:8)).*(grad_N(:,9:10)*TK),2));
    Aloc(4,6) = sum(QuadRule.w.*sum((grad_N(:,7:8)).*(grad_N(:,11:12)*TK),2));
  
    Aloc(5,5) = sum(QuadRule.w.*sum((grad_N(:,9:10)).*(grad_N(:,9:10)*TK),2));
    Aloc(5,6) = sum(QuadRule.w.*sum((grad_N(:,9:10)).*(grad_N(:,11:12)*TK),2));
  
    Aloc(6,6) = sum(QuadRule.w.*sum((grad_N(:,11:12)).*(grad_N(:,11:12)*TK),2));
   
  end

  % Fill in lower triangular part
  
  tri = triu(Aloc);
  Aloc = tri+tril(tri',-1);
 
return