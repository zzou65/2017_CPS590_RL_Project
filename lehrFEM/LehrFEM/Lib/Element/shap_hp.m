function shap = shap_hp(x,p)
% SHAP_HP Hierarchical shape functions.
%
%   SHAP = SHAP_HP(X,P) computes the values and gradients of the
%   hierarchical shape functions on the reference element up to polynomial
%   degree P at the points X.
%
%   Example:
%
%   Shap = shap_hp([0 0],10);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(x,1);
  aux = ones(1,2);
  
  % Preallocate memory
  
  shap.vshap.shap = cell(1,3);
  shap.vshap.grad_shap = cell(1,3);
  
  shap.eshap = cell(1,3);
  for i = 1:3
    shap.eshap{i}.shap = cell(1,p-1);
    shap.eshap{i}.grad_shap = cell(1,p-1);
  end
  
  shap.cshap.shap = cell(1,1);
  shap.cshap.grad_shap = cell(1,1);
  
  % Compute barycentric coordinates
  
  lam_1 = 1-x(:,1)-x(:,2);
  lam_2 = x(:,1);
  lam_3 = x(:,2);
  
  grad_lam_1 = -ones(nPts,2);
  grad_lam_2 = [ones(nPts,1) zeros(nPts,1)];
  grad_lam_3 = [zeros(nPts,1) ones(nPts,1)];
  
  % Compute vertex shape functions

  shap.vshap.shap{1} = lam_1;
  shap.vshap.grad_shap{1} = grad_lam_1;
  
  shap.vshap.shap{2} = lam_2;
  shap.vshap.grad_shap{2} = grad_lam_2;
  
  shap.vshap.shap{3} = lam_3;
  shap.vshap.grad_shap{3} = grad_lam_3;
  
  % Initialize auxiliary variables for edge shape functions
  
  y1 = lam_2-lam_3;
  grad_y1 = grad_lam_2-grad_lam_3;
  z1 = lam_2+lam_3;
  grad_z1 = grad_lam_2+grad_lam_3;
  v1 = ones(nPts,1);
  grad_v1 = zeros(nPts,2);
  w1 = y1;
  grad_w1 = grad_y1;
  
  y2 = lam_3-lam_1;
  grad_y2 = grad_lam_3-grad_lam_1;
  z2 = lam_3+lam_1;
  grad_z2 = grad_lam_3+grad_lam_1;
  v2 = ones(nPts,1);
  grad_v2 = zeros(nPts,2);
  w2 = y2;
  grad_w2 = grad_y2;
  
  y3 = lam_1-lam_2;
  grad_y3 = grad_lam_1-grad_lam_2;
  z3 = lam_1+lam_2;
  grad_z3 = grad_lam_1+grad_lam_2;
  v3 = ones(nPts,1);
  grad_v3 = zeros(nPts,2);
  w3 = y3;
  grad_w3 = grad_y3;
    
  for i = 0:(p-2)
      
    % Compute edge shape functions on edge 1 (opposite vertex 1)
    
    tmp_1 = w1;
    tmp_2 = grad_w1;
    grad_w1 = ((2*i+3)*((w1*aux).*grad_y1+(y1*aux).*grad_w1) - ...
               (i+1)*((2*z1.*v1*aux).*grad_z1+(z1.^2*aux).*grad_v1))/(i+2);
    w1 = ((2*i+3)*y1.*w1-(i+1)*z1.^2.*v1)/(i+2);
    shap.eshap{1}.shap{i+1} = 1/(2*i+3)*(w1-z1.^2.*v1);
    shap.eshap{1}.grad_shap{i+1} = 1/(2*i+3)*(grad_w1 - ...
                                              (2*z1.*v1*aux).*grad_z1 - ...
                                              (z1.^2*aux).*grad_v1);         
    v1 = tmp_1;
    grad_v1 = tmp_2;

    % Compute edge shape functions on edge 2 (opposite vertex 2)
        
    tmp_1 = w2;
    tmp_2 = grad_w2;
    grad_w2 = ((2*i+3)*((w2*aux).*grad_y2+(y2*aux).*grad_w2) - ...
               (i+1)*((2*z2.*v2*aux).*grad_z2+(z2.^2*aux).*grad_v2))/(i+2);
    w2 = ((2*i+3)*y2.*w2-(i+1)*z2.^2.*v2)/(i+2);
    shap.eshap{2}.shap{i+1} = 1/(2*i+3)*(w2-z2.^2.*v2);
    shap.eshap{2}.grad_shap{i+1} = 1/(2*i+3)*(grad_w2 - ...
                                              (2*z2.*v2*aux).*grad_z2 - ...
                                              (z2.^2*aux).*grad_v2);
    v2 = tmp_1;
    grad_v2 = tmp_2;
   
    % Compute edge shape functions on edge 3 (opposite vertex 3)
    
    tmp_1 = w3;
    tmp_2 = grad_w3;
    grad_w3 = ((2*i+3)*((w3*aux).*grad_y3+(y3*aux).*grad_w3) - ...
               (i+1)*((2*z3.*v3*aux).*grad_z3+(z3.^2*aux).*grad_v3))/(i+2);
    w3 = ((2*i+3)*y3.*w3-(i+1)*z3.^2.*v3)/(i+2);    
    shap.eshap{3}.shap{i+1} = 1/(2*i+3)*(w3-z3.^2.*v3);
    shap.eshap{3}.grad_shap{i+1} = 1/(2*i+3)*(grad_w3 - ...
                                              (2*z3.*v3*aux).*grad_z3 - ...
                                              (z3.^2*aux).*grad_v3);         
    v3 = tmp_1;
    grad_v3 = tmp_2;

  end
  
  % Initialize auxiliary variables for element shape functions
  
  v = zeros(nPts,p);
  grad_v = zeros(nPts,2*p);
 
  v(:,1) = lam_3;
  v(:,2) = lam_3.*(2*lam_3-1);
  grad_v(:,1:2) = grad_lam_3;
  grad_v(:,3:4) = (4*lam_3-1)*aux.*grad_lam_3;
  
  id = 1:2; 
  y = 2*lam_3-1;
  grad_y = 2*grad_lam_3;
  for i = 0:(p-1)
    v(:,i+3) = ((2*i+3)*y.*v(:,i+2)-(i+1)*v(:,i+1))/(i+2);
    grad_v(:,id+4) = ((2*i+3)*(v(:,i+2)*aux.*grad_y + ...
                               y*aux.*grad_v(:,id+2)) - ...
                      (i+1)*grad_v(:,id))/(i+2);
    id = id+2;
  end
  
  % Compute element shape functions
  
  k = 1;
  for i = 0:(p-3)
    for j = 0:i
      shap.cshap.shap{k} = shap.eshap{3}.shap{j+1}.*v(:,i-j+1);     
      shap.cshap.grad_shap{k} = shap.eshap{3}.grad_shap{j+1}.*(v(:,i-j+1)*aux) + ...
                                (shap.eshap{3}.shap{j+1}*aux).*grad_v(:,2*(i-j)+(1:2));
      k = k+1;
    end
  end
    
return