function Lieloc = LIEUP_FD(Vertices, dummy, vHandle, varargin)
% MASS_LFE Element mass matrix.
%
%   LIELOC = LIEUP_FD(VERTICES,v) computes the element Lie-Derivative matrix using linear
%   Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   dummy: useless variable, neccessary for interface to assemMat_LFE 
%
%  Vhandle handles velocity field
%
%
%   Example:
%
%   lieloc = LIEUP_FD(Vertices, v);

%   Copyright 2005-2005 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  Lieloc = zeros(3,3);
  B=zeros(3,3);
 
  a1 = Vertices(1,:);
  a2 = Vertices(2,:);
  a3 = Vertices(3,:);

  % Compute element mapping

  bK = a1;
  BK = [a2-bK; ...
        a3-bK];
  det_BK = abs(det(BK));
  inv_BK = inv(BK); 
  
  % Compute local mass matrix
  
  Mloc = det_BK/8*[2 0 0; 0 1 0; 0 0 1];
  
  % Compute gradient of local shape functions
  
  gradN1=[-1,-1]*inv_BK';
  gradN2=[1,0]*inv_BK';
  gradN3=[0,1]*inv_BK';
  
  % Extract nodal vectors
    
  v1 = -vHandle(a1);
  v2 = -vHandle(a2);
  v3 = -vHandle(a3);
  
  % Compute Lie derivative scheme for first vertex
  y = a1 + v1;
  yhat = (y-bK)*inv_BK;
    
  z = sum(yhat);
  if(z > 0)
   xhat = yhat/z;
    if(xhat(1) >= 0 && xhat(2) >= 0)
     B(1,:) = [(gradN2+gradN3)*v1', -gradN2*v1', -gradN3*v1'];
     if (xhat(1)==0 || xhat(2)==0)
      B(1,:)=1/2 * B(1,:);
     end
    end
  end    
           
  % Compute Lie derivative scheme for second vertex
      
  y = a2 + v2;
  yhat = (y-bK)*inv_BK;
    
  z = 1-yhat(1);
  if(z > 0)
   xhat = [0 yhat(2)/z];
    if(xhat(2) >= 0 && xhat(2) <= 1)
     B(2,:) = [-gradN1*v2', (gradN1+gradN3)*v2', -gradN3*v2'];
     if (xhat(2)==0 || xhat(1)== -xhat(2))
      B(2,:)=1/2 * B(2,:);
     end  
    end
  end
    
  % Compute Lie derivative scheme for third index
    
  y = a3 + v3;
  yhat = (y-bK)*inv_BK;
  
  z = 1-yhat(2);
  if(z > 0)
   xhat = [yhat(1)/z 0];
   if(xhat(1) >= 0 && xhat(1) <= 1)
    B(3,:) = [-gradN1*v3', -gradN2*v3', (gradN1+gradN2)*v3'];
    if (xhat(1)==0 || xhat(1) == -xhat(2))
     B(3,:)=1/2* B(3,:);
    end
   end
  end
 Lieloc=B;
return