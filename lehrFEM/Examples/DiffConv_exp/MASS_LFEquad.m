function Aloc =MASS_LFEquad(Vertices, flag, QuadRule)
% MASS_LFEquad Element mass matrix for the linear finite elements.
%
%   MLOC = MASS_LFEQUAD(VERTICES, flag, QuadRule) computes the element 
%   mass matrix using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Example:
%
%   Mloc = MASS_LFE(Vertices,1,P103());

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
%

  % Preallocate memory
  
  Aloc = zeros(3,3);

  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ];          % transpose of transformation matrix
  det_BK = abs(det(BK));               % twice the area of the triagle
  inv_BK=inv(BK);
  inv_BK_t=transpose(inv_BK);  
  
  % Gradient of Shapfunctions
  lambda = shap_LFE(QuadRule.x);
  
  w = QuadRule.w;
  
  for i=1:3
    for j=1:i
        Aloc(i,j)=sum(w.*lambda(:,i).*lambda(:,j))*det_BK;
    end
  end
Aloc=Aloc+tril(Aloc,-1)';
return