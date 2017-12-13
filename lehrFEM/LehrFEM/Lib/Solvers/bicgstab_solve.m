function varargout = bicgstab_solve(x,A,b,tol,maxit,varargin)
% BICGSTAB_SOLVE BiConjugate Gradient Stabilized Method. 
%
%   X = BICGSTAB_SOLVE(X,A,B,TOL,MAXIT) solves the linear system A*X=B from
%   the initial guess X0 using the BiCG stabilized method. 
%
%   The tolerance TOL and the maximum number of iterations MAXIT are used
%   as stopping criteria for the stabilized BiCG iterations.
%
%   [X,FLAG] = BICGSTAB_SOLVE(X0,A,B,TOL,MAXIT) the integer FLAG specifies
%   wheter the method has converged or not:
%    1 Method has converged within MAXIT steps to the prescribed tolerance
%      TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%
%   [X,FLAG,RELRES] = BICGSTAB_SOLVE(X0,A,B,TOL,MAXIT) also returns the
%   relative residual RELRES at the last stabilized BiCG iteration.
% 
%   [X,FLAG,RELRES,ITER] = BICGSTAB_SOLVE(X0,A,B,TOL,MAXIT) also returns
%   the number of stabilized BiCG iterations that were performed.
%
%   Example:
%
%   x = bicgstab_solve(zeros(size(b)),A,b,10,1e-6,1000);

%   Iterative template routine --
%   Univ. of Tennessee and Oak Ridge National Laboratory October 1, 1993
%   Details of this algorithm are described in
%   "Templates for the Solution of Linear Systems: Building Blocks for
%   Iterative Methods",
%   Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
%   and van der Vorst,
%   SIAM Publications, 1993.
%   (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).

  bnrm = norm(b);
  if(bnrm == 0.0)
    bnrm = 1.0;
  end
  
  % Run BiCG stabilized method
  
  if(isa(A,'function_handle'))
      
    % Stiffness matrix is given in terms of a function handle  
      
    r = b-A(x,varargin{:});
    rhat = r;
    alpha = 1.0;
    iter = 0;
    relres = tol+1;
    resvec = zeros(1,maxit);
    while(relres > tol & iter < maxit)
       rho_new = rhat'*r;
       if(rho_new == 0.0)
         error('BiCG stabilized method failed');
       end
       if(iter == 0)
         p = r;
       else
         beta  = (rho_new/rho_old)*(alpha/omega);
         p = r + beta*(p-omega*v);
       end
       v = A(p,varargin{:});
       alpha = rho_new/(rhat'*v);
       s = r - alpha*v;
       if(norm(s) < tol)
         x = x + alpha*p;
         iter = iter+1;
         relres = norm(s)/bnrm;
         resvec(iter) = relres;
         break;
       end
       t = A*s;
       omega = (t'*s)/(t'*t);
       if(omega == 0.0)
         error('BiCG stabilized method failed');
       end
       x = x + alpha*p + omega*s; 
       r = s - omega*t;
       rho_old = rho_new;
       iter = iter+1;
       relres = norm(r)/bnrm;
       resvec(iter) = relres;
    end

  else
    
    % Stiffness matrix is given in terms of a double matrix
    
    r = b-A*x;
    rhat = r;
    iter = 0;
    relres = tol+1;
    resvec = zeros(1,maxit);
    while(relres > tol & iter < maxit)
       rho_new = rhat'*r;
       if(rho_new == 0.0)
         error('BiCG stabilized method failed');
       end
       if(iter == 0)
         p = r;
       else
         beta  = (rho_new/rho_old)*(alpha/omega);
         p = r + beta*(p-omega*v);
       end
       v = A*p;
       alpha = rho_new/(rhat'*v);
       s = r - alpha*v;
       if(norm(s) < tol)
         x = x + alpha*p;
         iter = iter+1;
         relres = norm(s)/bnrm;
         resvec(iter) = relres;
         break;
       end
       t = A*s;
       omega = (t'*s)/(t'*t);
       if(omega == 0.0)
         error('BiCG stabilized method failed');
       end
       x = x + alpha*p + omega*s; 
       r = s - omega*t;
       rho_old = rho_new;
       iter = iter+1;
       relres = norm(r)/bnrm;
       resvec(iter) = relres;
    end
      
  end
  
  % Check convergence
  
  if(relres < tol)
    flag =  1;
  else
    flag = 0; 
  end

  % Assign output arguments
  
  if(nargout == 1)
    varargout{1} = x;  
  elseif(nargout == 2)
    varargout{1} = x;
    varargout{2} = flag;
  elseif(nargout == 3)
    varargout{1} = x;
    varargout{2} = flag;
    varargout{3} = relres;
  elseif(nargout == 4)
    varargout{1} = x;
    varargout{2} = flag;
    varargout{3} = relres;
    varargout{4} = iter;
  elseif(nargout == 5)
    varargout{1} = x;
    varargout{2} = flag;
    varargout{3} = relres;
    varargout{4} = iter;
    varargout{5} = resvec(1:iter);
  end
  
return
