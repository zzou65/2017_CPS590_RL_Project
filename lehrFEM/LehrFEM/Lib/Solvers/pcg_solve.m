function varargout = pcg_solve(x,A,b,tol,maxit,prec,varargin)
% PCG_SOLVE Preconditioned conjugate gradients method.
%
%   X = PCG_SOLVE(X0,A,B,TOL,MAXIT,PREC) computes the solution X to the
%   linear system A*X = B from the initial guess X0 using the CG method
%   together with the preconditioner PREC.
%
%   The tolerance TOL serves as a stopping criterion for the number of CG
%   iterations. MAXIT specifies the maximum number of CG iterations to be
%   used.
%
%   [X,FLAG] = PCG_SOLVE(X0,A,B,TOL,MAXIT,PREC) the integer FLAG specifies
%   wheter the method has converged or not:
%    1 Method has converged within MAXIT steps to the prescribed tolerance
%      TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%
%   [X,FLAG,RELRES] = PCG_SOLVE(X0,A,B,TOL,MAXIT,PREC) also returns the
%   relative residual RELRES at the last CG iteration.
% 
%   [X,FLAG,RELRES,ITER] = PCG_SOLVE(X0,A,B,TOL,MAXIT,PREC) also returns
%   the number of CG iterations that were performed.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = PCG_SOLVE(X0,A,B,TOL,MAXIT,PREC) also
%   returns the values of the relative residuals RESVEC at each CG
%   iteration.
%
%   Example:
%
%   x = pcg_solve(x0,A,b,1e-4,100,@bpx_prec);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  
      
  bnrm = norm(b);
  if(bnrm == 0.0)
    bnrm = 1;
  end

  % Run preconditioned CG solver
    
  if(isa(A,'function_handle'))
  
    % Stiffness matrix is given in terms of a function handle  
      
    r = b-A(x,varargin{:});
    iter = 0;
    relres = tol+1;
    resvec = zeros(1,maxit);
    while(iter < maxit & relres > tol)
      z = prec(r,varargin{:});
      rho_new = r'*z;
      if(iter == 0)
        p = z;  
      else
        beta = rho_new/rho_old;
        p = z + beta*p;
      end
      q = A(p,varargin{:});
      alpha = rho_new/(p'*q);
      x = x + alpha*p;
      r = r - alpha*q;
      rho_old = rho_new;
      iter = iter+1;
      relres = norm(r)/bnrm;
      resvec(iter) = relres;
    end
 
  else
     
    % Stiffness matrix is given in terms of double valued matrix  
      
    r = b-A*x;
    iter = 0;
    relres = tol+1;
    resvec = zeros(1,maxit);
    while(iter < maxit & relres > tol)
      z = prec(r,varargin{:});
      rho_new = r'*z;
      if(iter == 0)
        p = z;  
      else
        beta = rho_new/rho_old;
        p = z + beta*p;
      end
      q = A*p;
      alpha = rho_new/(p'*q);
      x = x + alpha*p;
      r = r - alpha*q;
      rho_old = rho_new;
      iter = iter+1;
      relres = norm(r)/bnrm;
      resvec(iter) = relres;
    end   
      
  end
  
  % Check convergence
  
  if(relres > tol)
    flag = 0;
  else
    flag = 1;  
  end
  
  % Assign output arguments
  
  switch(nargout)
    case 1
      varargout{1} = x;
    case 2
      varargout{1} = x;
      varargout{2} = flag;
    case 3
      varargout{1} = x;
      varargout{2} = flag;
      varargout{3} = relres;
    case 4
      varargout{1} = x;
      varargout{2} = flag;
      varargout{3} = relres;
      varargout{4} = iter;
    case 5
      varargout{1} = x;
      varargout{2} = flag;
      varargout{3} = relres;
      varargout{4} = iter;
      varargout{5} = resvec(1:iter);
  end
  
return