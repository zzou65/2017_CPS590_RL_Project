function Lieloc = STIMA_ContrGrad_Up(Vertices, dummy, vHandle, varargin)
% STIMA_ContrGrad_UP  Element matrix for convection term using upwind
% quadrature rule 
%
%   LIELOC = STIMA_ContrGrad_Up(VERTICES,d,v) computes the convection element 
%   matrix using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   dummy: useless variable, neccessary for interface to assemMat_LFE 
%
%   Vhandle: handles velocity field
%
%
%   Example:
%
%   lieloc = STIMA_ContrGrad_Up(Vertices,0,v);
%
%   Copyright 2007 Patrick Meury, Holger Heumann
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
  
  % Compute gradient of local shape functions
  
  gradN1=[-1,-1]*inv_BK';
  gradN2=[1,0]*inv_BK';
  gradN3=[0,1]*inv_BK';
  
  % Extract nodal vectors
    
  v1 = -vHandle(a1,varargin{:});
  v2 = -vHandle(a2,varargin{:});
  v3 = -vHandle(a3,varargin{:});
  
  % Compute Lie derivative scheme for first vertex using upwind quadrature
  % rule
  y = a1 + v1;
  yhat = (y-bK)*inv_BK;
    
  z = sum(yhat);
  if(z > 0)
   xhat = yhat/z;
    if(xhat(1) >= 0 && xhat(2) >= 0)
     B(1,:) = [(gradN2+gradN3)*v1', -gradN2*v1', -gradN3*v1'];
     if (xhat(1)==0 || xhat(2)==0)   % two elements contribute to matrix
      B(1,:)=1/2 * B(1,:);
     end
    end
  end    
           
  % Compute Lie derivative scheme for second vertex using upwind quadrature
  % rule
      
  y = a2 + v2;
  yhat = (y-bK)*inv_BK;
    
  z = 1-yhat(1);
  if(z > 0)
   xhat = [0 yhat(2)/z];
    if(xhat(2) >= 0 && xhat(2) <= 1)
     B(2,:) = [-gradN1*v2', (gradN1+gradN3)*v2', -gradN3*v2'];
     if (xhat(2)==0 || xhat(1)== -xhat(2)) 
      B(2,:)=1/2 * B(2,:);   %two elements contribute to matrix
     end  
    end
  end
    
  % Compute Lie derivative scheme for third index using upwind quadrature
    
  y = a3 + v3;
  yhat = (y-bK)*inv_BK;
  
  z = 1-yhat(2);
  if(z > 0)
   xhat = [yhat(1)/z 0];
   if(xhat(1) >= 0 && xhat(1) <= 1)
    B(3,:) = [-gradN1*v3', -gradN2*v3', (gradN1+gradN2)*v3'];
    if (xhat(1)==0 || xhat(1) == -xhat(2))
     B(3,:)=1/2* B(3,:);    %two elements contribute to matrix 
    end
   end
  end
 Lieloc=B;
return