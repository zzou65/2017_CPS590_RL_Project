function [MD_loc, UP_loc, C_loc] = STIMA_UPQFE(Vertices,vHandle, weights, varargin)
% STIMA_UPQFE  Element matrix for convection term using upwind
% quadrature rule with weights at vertices, midpoints and barycenter
%
%   LIELOC = STIMA_UPQFE(VERTICES,v) computes the convection element 
%   matrix using bilinear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   Vhandle: handles velocity field
%    
%    weights=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45]
%    weights=[0 0 0 1/3 1/3 1/3 0];
%    weights=[1/12 1/12 1/12 0 0 0 3/4];
%
%   Example:
%
%   lieloc = STIMA_UPQFE(Vertices,0,v);
%
%   Copyright 2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  UP_loc = zeros(6,6);
  C_loc=zeros(6,6);
 
  % Quadrature weights and points
  %w_c=0.45  barycenter; w_v=0.05  vertices
  %w_m=4/30  midpoints on edges
  
%  w=[0.05 0.05 0.05 4/30 4/30 4/30];
  w=weights(1:6);
  %w=[0 0 0 1/3 1/3 1/3];

  x=[ 0 0; 1 0; 0 1;0.5 0; 0.5 0.5; 0 0.5];
  %w_c=0.45;
  w_c=weights(7);
  x_c=[1/3 1/3];
  % Compute element mapping

  a1 = Vertices(1,:);
  a2 = Vertices(2,:);
  a3 = Vertices(3,:);
  
  m1=(a2+a1)/2;
  m2=(a3+a2)/2;
  m3=(a1+a3)/2;
  
  b=(a1+a2+a3)/3;
  
  bK = a1;
  BK = [a2-bK; ...
        a3-bK];
  det_BK = abs(det(BK));
  inv_BK = inv(BK); 
  
  % Compute gradient of local shape functions
  
  gradN=grad_shap_QFE(x);
  for i=1:2:12
      gradN(:,[i i+1])=gradN(:,[i i+1])*inv_BK';
  end
  
  % Extract nodal vectors
    
  v = -vHandle(x*BK + ones(6,1)*a1);
  
  % Compute Lie derivative scheme for first vertex using upwind quadrature
  % rule
  y = a1 + v(1,:);
  yhat = (y-bK)*inv_BK;
    
  z = sum(yhat);
  if(z > 0)
   xhat = yhat/z;
    if(xhat(1) >= 0 && xhat(2) >= 0)
     UP_loc(1,:) = [-gradN(1,[1 2])*v(1,:)', -gradN(1,[3 4])*v(1,:)', ....
                        -gradN(1,[5 6])*v(1,:)', -gradN(1,[7 8])*v(1,:)', ....
                        -gradN(1,[9 10])*v(1,:)', -gradN(1,[11 12])*v(1,:)'];
     if (xhat(1)==0 || xhat(2)==0)   % two elements contribute to matrix
      UP_loc(1,:)=1/2 * UP_loc(1,:);
     end
    end
  end    
           
  % Compute Lie derivative scheme for second vertex using upwind quadrature
  % rule
      
  y = a2 + v(2,:);
  yhat = (y-bK)*inv_BK;
    
  z = 1-yhat(1);
  if(z > 0)
   xhat = [0 yhat(2)/z];
    if(xhat(2) >= 0 && xhat(2) <= 1)
     UP_loc(2,:) =[-gradN(2,[1 2])*v(2,:)', -gradN(2,[3 4])*v(2,:)', ...
                        -gradN(2,[5 6])*v(2,:)', -gradN(2,[7 8])*v(2,:)', ...
                         -gradN(2,[9 10])*v(2,:)', -gradN(2,[11 12])*v(2,:)'];
     if (xhat(2)==0 || xhat(1)== -xhat(2)) 
      UP_loc(2,:)=1/2 * UP_loc(2,:);   %two elements contribute to matrix
     end  
    end
  end
    
  % Compute Lie derivative scheme for third index using upwind quadrature
    
  y = a3 + v(3,:);
  yhat = (y-bK)*inv_BK;
  
  z = 1-yhat(2);
  if(z > 0)
   xhat = [yhat(1)/z 0];
   if(xhat(1) >= 0 && xhat(1) <= 1)
    UP_loc(3,:) =[-gradN(3,[1 2])*v(3,:)', -gradN(3,[3 4])*v(3,:)', ...
                               -gradN(3,[5 6])*v(3,:)', -gradN(3,[7 8])*v(3,:)', ...
                               -gradN(3,[9 10])*v(3,:)', -gradN(3,[11 12])*v(3,:)'];
    if (xhat(1)==0 || xhat(1) == -xhat(2))
     UP_loc(3,:)=1/2* UP_loc(3,:);    %two elements contribute to matrix 
    end
   end
  end

  %Compute Lie derivative scheme for first midpoint using upwind quadrature    
  yhat = (m1+v(4,:)-bK)*inv_BK;
  if(yhat(2) >= 0)
    UP_loc(4,:) = [-gradN(4,[1 2])*v(4,:)',   -gradN(4,[3 4])*v(4,:)', ...
                        -gradN(4,[5 6])*v(4,:)',   -gradN(4,[7 8])*v(4,:)', ...
                        -gradN(4,[9 10])*v(4,:)', -gradN(4,[11 12])*v(4,:)'];
  end
  if (yhat(2)==0)
     UP_loc(4,:)=1/2* UP_loc(4,:);    %two elements contribute to matrix 
  end
  
  % Compute Lie d6erivative scheme for second midpoint using upwind quadrature    
  yhat = (m2+v(5,:)-bK)*inv_BK;
  if(sum(yhat) <= 1)
    UP_loc(5,:) = [-gradN(5,[1 2])*v(5,:)', -gradN(5,[3 4])*v(5,:)', ...
                       -gradN(5,[5 6])*v(5,:)', -gradN(5,[7 8])*v(5,:)', ...
                       -gradN(5,[9 10])*v(5,:)', -gradN(5,[11 12])*v(5,:)'];
  end
  if (sum(yhat)==1)
     UP_loc(5,:)=1/2* UP_loc(5,:);    %two elements contribute to matrix 
  end

  % Compute Lie derivative scheme for third midpoint using upwind quadrature
  yhat = (m3+v(6,:)-bK)*inv_BK;
  if(yhat(1) >= 0)
    UP_loc(6,:) = [-gradN(6,[1 2])*v(6,:)', -gradN(6,[3 4])*v(6,:)', ...
                        -gradN(6,[5 6])*v(6,:)', -gradN(6,[7 8])*v(6,:)', ...
                        -gradN(6,[9 10])*v(6,:)', -gradN(6,[11 12])*v(6,:)'];
  end
  if (yhat(1)==0)
     UP_loc(6,:)=1/2* UP_loc(6,:);    %two elements contribute to matrix 
  end

  % compute contribution of center quadrature point
  v_c=vHandle(b);
  N_c=shap_QFE(x_c);
  grad_N_c=grad_shap_QFE(x_c);
  grad_N_c=inv_BK*reshape(grad_N_c,2,6);
  v_gradN_c=sum(grad_N_c.*(v_c'*[1 1 1 1 1 1]),1);
  
  C_loc=det_BK/2*w_c*v_gradN_c'*N_c;
  %C_loc=det_BK/2*v_gradN_c'*N_c;
  C_loc=C_loc';
  
  % output
  MD_loc=det_BK/2*w;
return
