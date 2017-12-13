function x = curlauxpc(r,A,P,VL,G,Lpot,m)
% auxiliary space preconditioner for H(curl)-elliptic problem
% discretized by means of edge elements
%
% r -> (residual type) argument passed to the preconditioner
% A -> edge element Galerkin matrix
% P -> nodal-to-edge transfer matrix
% VL -> Galerkin matrix in nodal auxiliary space
% G -> discrete gradient matrix
% Lpot -> projection of A to discrete potential space
% m -> number of Gauss-Seidel smoothing steps (optional)

% DEBUGGING CODE
% fprintf('size(A) = '); disp(size(A));
% fprintf('size(P) = '); disp(size(P));
% fprintf('size(VL) = '); disp(size(VL));
% fprintf('size(G) =  '); disp(size(G));
% fprintf('size(Lpot) = '); disp(size(Lpot));

if (nargin < 7), m = 1; end

x = zeros(length(r),1);

if (m > 0)
% m steps of symmetric Gauss-Seidel
% if m=0, then use Jacobi smoother

  Le = tril(A);
  Re = triu(A);
  
  for i = 1:m, x = x + Le\(r-A*x); end
  for i = 1:m, x = x + Re\(r-A*x); end
else
  x = r./diag(A);
end

% Additive auxiliary space correction on Lagrangian finite element space
x = x + P*(VL\(transpose(P)*r));

% Additive auxiliary space correction in discrete potential space
x = x +  G*(Lpot\(G'*r));

end
