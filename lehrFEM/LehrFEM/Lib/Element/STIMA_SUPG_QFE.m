function Aloc = STIMA_SUPG_QFE(Vertices, flag, QuadRule, VHandle, a,d1,d2, varargin)
% STIMA_SUPG_QFE Element stiffness matrix for the Laplacian.
%
%   ALOC = STIMA_SUPG_QFE(VERTICES) computes the element stiffness matrix
%   for SUPG-modification for convection using linear Lagrangian finite elements.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%  
%   a: diffusivity
%   d1 d2: apriori chosen constants for SUPG-modification
%
%   Flag useless, needed for interface to assemMat_LFE
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.e: 
%
%   VHANDLE is function handle for velocity field   
%
%   Example:
%
%   Aloc = STIMA_SUPG_QFE([0 0; 1 0; 0 1]);

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
  
  % Quadrature points in actual element stored as rows of a matrix
  nPoints = size(QuadRule.w,1);
  x = QuadRule.x*BK + ones(nPoints,1)*P1;
  
  % first order penalty
  % Gradient of Shapfunctions
  grad_lambda = grad_shap_QFE(QuadRule.x);
  
  for k=1:2:12
      grad_lambda(:,k:k+1)=grad_lambda(:,k:k+1)*inv_BK_t;
  end

  % Evaluate coefficient function at quadrature nodes
  c =VHandle(x,varargin{:});
  
  % velocity times grad basisfunctions
  
  w = QuadRule.w;
  cgl=zeros(nPoints,6);
  for k=1:2:12
    cgl(:,(k+1)/2)=sum(c.*grad_lambda(:,k:k+1),2);
  end
  
for i=1:6
    for j=1:i
        Aloc(i,j)=sum(w.*cgl(:,i).*cgl(:,j))*det_BK;
    end
end
Aloc=Aloc+tril(Aloc,-1)';
  
  
  % second order penalty
  % Lapl of Shapfunctions
  ll_m=zeros(1,6);
  gl = grad_shap_LFE([1/3,1/3]); % constants
  gl(:,1:2)=gl(:,1:2)*inv_BK_t;
  gl(:,3:4)=gl(:,3:4)*inv_BK_t;
  gl(:,5:6)=gl(:,5:6)*inv_BK_t;
  ll_m(1)=4*gl(:,1:2)*gl(:,1:2)';
  ll_m(2)=4*gl(:,3:4)*gl(:,3:4)';
  ll_m(3)=4*gl(:,5:6)*gl(:,5:6)';
  ll_m(4)=8*gl(:,1:2)*gl(:,3:4)';
  ll_m(5)=8*gl(:,3:4)*gl(:,5:6)';
  ll_m(6)=8*gl(:,5:6)*gl(:,1:2)';
  
  aloc_l=zeros(1,6);
  for i=1:6
      aloc_l(i)=det_BK*sum(w.*cgl(:,i));
  end
  Aloc=Aloc-a*aloc_l'*ll_m;

  % local PecletNumber
  hK=max([norm(P2-P1),norm(P3-P1),norm(P2-P3)]);
  v_infK=max(abs(c(:)));
  PK=v_infK*hK/(2*a);
  
  if (PK<=1)
    Aloc=d1*hK^2/a*Aloc;
  else
    Aloc=d2*hK*Aloc;
  end
return