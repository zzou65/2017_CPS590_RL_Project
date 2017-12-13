function Aloc=STIMA_ConvInt_LFE(Vertices,flag, VHandle,QuadRule)
% local stiffness matrix for convetion term using Edge-Element Interpolation
% for velocity

%  Vertices: coordinates of vertices or local element
%  QuadRule: 1-D quadrule
%  VHandle: function handle to velocity

  % Preallocate memory
  
  Aloc=zeros(3,3);
  grad_lambda=zeros(1,6); 
  w=zeros(3,1);  
  
  nGauss = size(QuadRule.w,1);
  % Compute element mapping
  
  P1 = Vertices(1,:);
  P2 = Vertices(2,:);
  P3 = Vertices(3,:);
  
  BK = [ P2 - P1 ; P3 - P1 ];          % transpose of transformation matrix
  det_BK = abs(det(BK));               % twice the area of the triagle
  inv_BK=inv(BK);
  inv_BK_t=transpose(inv_BK);
  
  % gradients of  LFE-elements
  grad_lambda=grad_shap_LFE([1/3 1/3]);
  grad_lambda=reshape(inv_BK*reshape(grad_lambda,2,3),1,6);
  
  % Edge-Element-Interpolation  of velocity
  
  % first edge
  x = ones(nGauss,1)*P2+QuadRule.x*(P3-P2);
  dS = ones(nGauss,1)*(P3-P2);
  Fval = VHandle(x);
  w(1) = sum(QuadRule.w.*sum(Fval.*dS,2));
  % second edge
  x = ones(nGauss,1)*P3+QuadRule.x*(P1-P3);
  dS = ones(nGauss,1)*(P1-P3);
  Fval = VHandle(x);
  w(2) = sum(QuadRule.w.*sum(Fval.*dS,2));
  % third edge
  x = ones(nGauss,1)*P1+QuadRule.x*(P2-P1);
  dS = ones(nGauss,1)*(P2-P1);
  Fval = VHandle(x);
  w(3) = sum(QuadRule.w.*sum(Fval.*dS,2));
  
  w=[w(2)-w(3) w(3)-w(1) w(1)-w(2)];
  
  for i=1:2:6
      for j=i:2:6
          Aloc((i+1)/2,(j+1)/2)=grad_lambda(i:i+1)*grad_lambda(j:j+1)';
      end
  end

  Aloc=Aloc+triu(Aloc,1)';
  
  w=w*Aloc;
   
  Aloc=det_BK/18*[1 1 1]'*w;
  
end
  