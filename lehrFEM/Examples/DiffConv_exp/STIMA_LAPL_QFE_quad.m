function Aloc = STIMA_LAPL_QFE_quad(Vertices, flag, QuadRule)
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
  
  % Gradient of Shapfunctions
  grad_lambda = grad_shap_QFE(QuadRule.x);
  
  for k=1:2:12
      grad_lambda(:,k:k+1)=grad_lambda(:,k:k+1)*inv_BK_t;
  end
  
  w = QuadRule.w;
  
  for i=1:2:12
    for j=1:2:i
        Aloc((i+1)/2,(j+1)/2)=sum(w.*sum(grad_lambda(:,i:i+1).*grad_lambda(:,j:j+1),2))*det_BK;
    end
  end
Aloc=Aloc+tril(Aloc,-1)';
return