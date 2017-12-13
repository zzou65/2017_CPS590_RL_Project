function varargout = eigen_Lanczos(A,tol,maxit,varargin)
% EIGEN_LANCZOS Estimates the eigen value.
%
%   [MIN,MAX] = EIGEN_LANCZOS(A,TOL,MAXIT) approximates the largest and 
%   smallest eigenvalue within MAXIT steps up to the prescribed tolerance 
%   TOL using a Lanczos process.
%
%   [MIN,MAX] = EIGEN_LANCZOS(A,TOL,MAXIT,P) computes the largest and 
%   smallest eigenvalues of the preconditioned matrix P*A. 
%
%   [MIN,MAX] = EIGEN_LANCZOS(A,TOL,MAXIT,P,PARAM) also handles the variable
%   length argument list PARAM to the preconditioner P of the matrix A.
%
%   [MIN,MAX,FLAG] = EIGEN_LANCZOS(A,TOL,MAXIT) also returns the convergence
%   flag FLAG
%    0 Method converged within MAXIT steps to prescribed tolerance TOL.
%    1 Method failed to converge.
%
%   [MIN,MAX,FLAG,ITER] = EIGEN_LANCZOS(A,TOL,MAXIT) also return the
%   iteration number ITER at which the method converged.

  % Set up preconditioner
  
  if(nargin > 3)
    P = varargin{1};
    Param = varargin(2:end);
  else
    P = @(x)x;
    Param = {};
  end

  % Minimal number of iterations
  minit = 10;
  
  % Initialize
  
  w = rand(size(A,1),1);
  v = A*w;
  norm = sqrt(v'*w);
  v = P(v/norm,Param{:});
  w = w/norm;
  alpha(1) = w'*A*v;
  lam_max = alpha(1);
  lam_min = alpha(1);
  lam_max_old = 2*lam_max;
  lam_min_old = 2*lam_min;

  % Start Lanczos iteration
  
  k = 1;
  delta = inf;
  while((k < minit) || (k < maxit && delta > tol))
    v = v-alpha(k)*w;
    beta(k) = sqrt(v'*A*v);
    if(beta(k) < tol)
      break  
    end
    tmp = w;
    w = v/beta(k);
    v = -beta(k)*tmp;
    v = v+P(A*w,Param{:});
    k = k+1;
    alpha(k) = w'*A*v;
    beta(k) = 0;
    lam_min_old = lam_min;
    lam_max_old = lam_max;
    [lam_min,lam_max] = bisection(lam_min,lam_max,alpha,beta,tol);
    delta = max((lam_max-lam_max_old)/lam_max_old, ...
                (lam_min_old-lam_min)/lam_min);
  end
  
  % Check convergence
  
  if(k >= maxit)
    flag = 0;  
  else
    flag = 1;  
  end
  
  % Assign output arguments
  
  switch(nargout)
    case 2
      varargout{1} = lam_min;
      varargout{2} = lam_max;
    case 3
      varargout{1} = lam_min;
      varargout{2} = lam_max;
      varargout{3} = flag;
    case 4
      varargout{1} = lam_min;
      varargout{2} = lam_max;
      varargout{3} = flag;
      varargout{4} = k;
    otherwise
      error('Not enough output arguments');   
  end
        
return
  
%%% Bisection aglorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lam_min,lam_max] = bisection(lam_min,lam_max,alpha,beta,tol)

% Start bisection algorithm

  y_max = max(alpha + abs(beta) + [0,abs(beta(1:size(beta,2)-1))]);
  y_min = min(alpha - abs(beta) - [0,abs(beta(1:size(beta,2)-1))]);

  % Determine upper bound

  pz = polynom(y_max,alpha,beta);
  while(abs(y_max-lam_max)/min(abs(y_max),abs(lam_max)) > tol)
    x = 0.5*(y_max+lam_max);
    px = polynom(x,alpha,beta);
    if (px*pz <= 0)
      lam_max = x;
    else
      y_max = x;
      pz = px;
    end
  end
  py = polynom(lam_max,alpha,beta);
  if (py*pz <= 0) & (py ~= 0) 
    lam_max = y_max;
  end

  % Determine lower bound

  py = polynom(y_min,alpha,beta);
  while(abs(y_min-lam_min)/min(abs(y_min),abs(lam_min)) > tol)
    x = 0.5*(y_min+lam_min);
    px = polynom(x,alpha,beta);
    if (px*py <= 0)
      lam_min = x;
    else
      y_min = x;
      py = px;
    end
  end
  pz = polynom(lam_min,alpha,beta);
  if (py*pz <= 0) & (pz ~= 0)
    lam_min = y_min;
  end

return

%%% Polynomial subroutine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = polynom(x,alpha,beta)

  r = 1;
  p = alpha(1)-x;
  for l=2:size(alpha,2)
    q = p;
    p = (alpha(l)-x)*p-beta(l-1)^2*r;
    r = q;
  end
 
return 