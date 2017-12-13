function Aloc = STIMA_Grad_MIXDG(Vertices,ElemInfo,varargin)
% STIMA_GRAD_MIXDG Element stiffness matrix for the gradient of multiplier p.
%
%   ALOC = STIMA_GRAD_MIXDG(VERTICES,ELEMINFO) computes the element stiffness 
%   matrix for the gradient of multiplier p coupling discontinuous and
%   conforming nodal elements
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  QuadRule=P3O3();
  Aloc = zeros(6,3);
  
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  bK = P1;
  BK = [P2-bK;P3-bK];
  det_BK = abs(det(BK));
  
  N = shap_DGLFE(QuadRule.x);
  
  % Compute the discrete gradient
  
  gradient(1,1:2) = [P2(2)-P3(2) -P2(1)+P3(1)];
  gradient(2,1:2) = [P3(2)-P1(2) -P3(1)+P1(1)];
  gradient(3,1:2) = [P1(2)-P2(2) -P1(1)+P2(1)];
  
  % Compute element stiffness matrix

for i = 1:3
    for j = 1:3
        Aloc(2*i-1,j) = sum(QuadRule.w.*N(:,i))*gradient(j,1);
        Aloc(2*i,j) = sum(QuadRule.w.*N(:,i))*gradient(j,2);
    end
end
        
return