  B=zeros(3,6);
  QuadRule=P3O3();
  s=size(QuadRule.x,1);
  n=shap_LFE(QuadRule.x);
  for j=1:s
         B=B+QuadRule.w(j)*n(j,:)'*[n(j,1)*(2*n(j,1)-1) n(j,2)*(2*n(j,2)-1) n(j,3)*(2*n(j,3)-1) 4*n(j,1)*n(j,2) 4*n(j,2)*n(j,3) 4*n(j,1)*n(j,3)];
  end
  B
  C=1/120*[2 -1 -1 8 4 8; ...
                            -1 2 -1 8 8 4; ...
                            -1 -1 2 4 8 8]