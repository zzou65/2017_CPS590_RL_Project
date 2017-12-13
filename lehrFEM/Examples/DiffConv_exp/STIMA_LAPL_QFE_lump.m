function Aloc = STIMA_LAPL_QFE_lump(Vertices, flag, varargin)
% STIMA_LAPL_QFE_lump Element stiffness matrix for the Laplacian.
% using barycentric lumping
%

  % Preallocate memory
  
  Aloc = zeros(6,6);

  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ];          % transpose of transformation matrix
  det_BK = abs(det(BK));               % twice the area of the triagle
  inv_BK=inv(BK);
  inv_BK_t=transpose(inv_BK);  
  
  %  Gradient of LFE-shape functions
  gl = grad_shap_LFE([1/3,1/3]);  % constants
  lambda1=gl(:,1:2)*inv_BK_t;
  lambda2=gl(:,3:4)*inv_BK_t;
  lambda3=gl(:,5:6)*inv_BK_t;
  
  % lumped element contribution
  
  Aloc(1,1)=9*lambda1*lambda1';
  Aloc(2,2)=9*lambda2*lambda2';
  Aloc(3,3)=9*lambda3*lambda3';
  Aloc(4,4)=16*(lambda1*lambda1'+lambda2*lambda2');
  Aloc(5,5)=16*(lambda2*lambda2'+lambda3*lambda3');
  Aloc(6,6)=16*(lambda1*lambda1'+lambda3*lambda3');
  
  Aloc(4,1)=12*lambda1*lambda2';
  Aloc(6,1)=12*lambda1*lambda3';
  Aloc(4,2)=12*lambda2*lambda1';
  Aloc(5,2)=12*lambda2*lambda3';
  Aloc(5,3)=12*lambda3*lambda2';
  Aloc(6,3)=12*lambda3*lambda1';
  
  Aloc(1:3,4:6)=Aloc(4:6,1:3)';
  
  Aloc(5,4)=8*lambda1*lambda3';
  Aloc(6,4)=8*lambda1*lambda2';
  Aloc(6,5)=8*lambda2*lambda3';

  Aloc(4,5)=Aloc(5,4);
  Aloc(4,6)=Aloc(6,4);
  Aloc(5,6)=Aloc(6,5);
  
  Aloc=det_BK/6*Aloc;
return