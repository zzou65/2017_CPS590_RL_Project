function UP_loc = STIMA_UPLFE(Vertices,vHandle, mass, varargin)
% STIMA_UPQFE  Element matrix for convection term using upwind
% quadrature rule 
%
%   LIELOC = STIMA_UPQFE(VERTICES,v, mass) computes the convection element 
%   matrix using bilinear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%    
%   MASS Mass of elements charing midpoints
%
%   Vhandle: handles velocity field
%    
%   Example:
%
%   lieloc = STIMA_UPQFE(Vertices,0,v);
%
%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  UP_loc = zeros(3,3);
 
  x=[0.5 0; 0.5 0.5; 0 0.5];
 
  % Compute element mapping

  a1 = Vertices(1,:);
  a2 = Vertices(2,:);
  a3 = Vertices(3,:);
  
  m1=(a2+a1)/2;
  m2=(a3+a2)/2;
  m3=(a1+a3)/2;
    
  bK = a1;
  BK = [a2-bK; ...
        a3-bK];
  det_BK = abs(det(BK));
  inv_BK = inv(BK); 
  
  % Compute gradient of local shape functions
  
  gradN=grad_shap_LFE(x);
  for i=1:2:6
      gradN(:,[i i+1])=gradN(:,[i i+1])*inv_BK';
  end
  
  % Extract nodal vectors
    
  v = -vHandle(x*BK + ones(3,1)*a1);
  
  %Compute Lie derivative scheme for first midpoint using upwind quadrature    
  yhat = (m1+v(1,:)-bK)*inv_BK;
  if(yhat(2) >= 0)
      elem= mass(1)/2*[-gradN(1,[1 2])*v(1,:)',...
                 -gradN(1,[3 4])*v(1,:)',...
                 -gradN(1,[5 6])*v(1,:)'];
    UP_loc(1,:) = elem; 
    UP_loc(2,:) = elem;
  end
  if (yhat(2)==0)
     UP_loc(1,:)=1/2* UP_loc(1,:);
     UP_loc(2,:)=1/2* UP_loc(1,:); %two elements contribute to matrix (what if edge on the boundary??????)
  end
  
  % Compute Lie derivative scheme for second midpoint using upwind quadrature    
  yhat = (m2+v(2,:)-bK)*inv_BK;
  if(sum(yhat) <= 1)
    elem = mass(2)/2*[-gradN(2,[1 2])*v(2,:)', ...
                -gradN(2,[3 4])*v(2,:)', ...
                -gradN(2,[5 6])*v(2,:)'];
    UP_loc(3,:) = UP_loc(3,:)+elem;
    UP_loc(2,:) = UP_loc(2,:)+elem;
  end
  if (sum(yhat)==1)
     UP_loc(2,:)=1/2* UP_loc(2,:);
     UP_loc(3,:)=1/2* UP_loc(3,:);%two elements contribute to matrix 
  end

  % Compute Lie derivative scheme for third midpoint using upwind quadrature
  yhat = (m3+v(3,:)-bK)*inv_BK;
  if(yhat(1) >= 0)
    elem = mass(3)/2*[-gradN(3,[1 2])*v(3,:)', ... 
                -gradN(3,[3 4])*v(3,:)', ...
                -gradN(3,[5 6])*v(3,:)'];
    UP_loc(3,:) = UP_loc(3,:)+elem;
    UP_loc(1,:) = UP_loc(1,:)+elem;
  end
  if (yhat(1)==0)
     UP_loc(1,:)=1/2* UP_loc(1,:);
     UP_loc(3,:)=1/2* UP_loc(3,:); %two elements contribute to matrix
  end

return