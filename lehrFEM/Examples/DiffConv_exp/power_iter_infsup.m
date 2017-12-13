function varargout=power_iter_infsup(A,V,U,tol,maxit);

% power_iter_infsup(A,V,U,tol) calculates minimal generalized
% eigenvalue and eigenvector: A' V^-1 A' x =l U x; for Sangallis
% infsup constant.

[m,n] = size(A);
[L,R] = lu(A);
z_old = rand(m,1);
z_old = z_old/norm(z_old);
for i = 1:maxit
   z_new = U*z_old;
   z_new = transpose(R)\z_new;
   z_new = transpose(L)\z_new;
   z_new = V*z_new;
   z_new = L\z_new;
   z_new = R\z_new;
   l=norm(z_new);
   z_new=z_new/norm(z_new);
   if (norm(z_new-z_old)<tol)
       break;
   end
   z_old=z_new;
end

if (nargout==1)
    varargout{1}=1/l;
end

if (nargout==2)
    varargout{1}=1/l;
    varargout{2}=z_new;
end

if (nargout==3)
    varargout{1}=1/l;
    varargout{2}=z_new;
    varargout{3}=i;
end

