function Aloc =MASS_QFEquad(Vertices, flag, QuadRule)
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
  lambda = shap_QFE(QuadRule.x);
  
  w = QuadRule.w;
  
  for i=1:6
    for j=1:i
        Aloc(i,j)=sum(w.*lambda(:,i).*lambda(:,j))*det_BK;
    end
  end
Aloc=Aloc+tril(Aloc,-1)';
return