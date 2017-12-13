function varargout = gmres_solve(x,A,b,restart,tol,maxit,varargin)
% GMRES_SOLVE Generalized minimal residuals method. 
%
%   X = GMRES_SOLVE(X,A,B,RESTART,TOL,MAXIT) solves the linear system A*X=B
%   from the initial guess X0 using the GMRES method. 
%
%   The tolerance TOL and the maximum number of iterations MAXIT are used
%   as stopping criteria for the GMRES iterations.
%
%   [X,FLAG] = GMRES_SOLVE(X0,A,B,RESTART,TOL,MAXIT) the integer FLAG
%   specifies wheter the method has converged or not:
%    1 Method has converged within MAXIT steps to the prescribed tolerance
%      TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%
%   [X,FLAG,RELRES] = GMRES_SOLVE(X0,A,B,RESTART,TOL,MAXIT) also returns
%   the relative residual RELRES at the last GMRES iteration.
% 
%   [X,FLAG,RELRES,ITER] = GMRES_SOLVE(X0,A,B,RESTART,TOL,MAXIT) also
%   returns the number of GMRES iterations that were performed.
%
%   Example:
%
%   x = gmres_solve(zeros(size(b)),A,b,10,1e-6,1000);

%   Iterative template routine --
%   Univ. of Tennessee and Oak Ridge National Laboratory October 1, 1993
%   Details of this algorithm are described in
%   "Templates for the Solution of Linear Systems: Building Blocks for
%   Iterative Methods",
%   Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
%   and van der Vorst,
%   SIAM Publications, 1993.
%   (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).

  % Initialize counters

  bnrm = norm(b);
  if(bnrm == 0.0)
    bnrm = 1.0;
  end
  n = size(b,1);
  m = restart;
  
  % Initialize workspace
  
  V = zeros(n,m+1);
  H = zeros(m+1,m);
  cs = zeros(m,1);
  sn = zeros(m,1);
  e1 = zeros(n,1);
  e1(1) = 1.0;

  % Run GMRES solver
  
  if(isa(A,'function_handle'))
  
    % Stiffness matrix is given in terms of a function handle  
      
    resvec = zeros(1,maxit);
    for iter = 1:maxit
      
      r = b-A(x,varargin{:});
      V(:,1) = r/norm(r);
      s = norm(r)*e1;
    
      % Construct orthonormal basis using Gram-Schmidt
    
      for i = 1:m
        w = A(V(:,i),varargin{:});
	    for k = 1:i,
	      H(k,i)= w'*V(:,k);
	      w = w-H(k,i)*V(:,k);
        end
	    H(i+1,i) = norm(w);
	    V(:,i+1) = w/H(i+1,i);
     
        % Apply Givens rotations
     
	    for k = 1:i-1
          temp = cs(k)*H(k,i)+sn(k)*H(k+1,i);
          H(k+1,i) = -sn(k)*H(k,i)+cs(k)*H(k+1,i);
          H(k,i) = temp;
        end
	    [cs(i),sn(i)] = rotmat(H(i,i),H(i+1,i)); 
      
        % Form i-th rotation matrix
      
        temp = cs(i)*s(i);                       
      
        % Approximate residual norm
      
        s(i+1) = -sn(i)*s(i);
	    s(i) = temp;
        H(i,i) = cs(i)*H(i,i)+sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
	    relres = abs(s(i+1))/bnrm;
     
        % Udpate approximation and exit
     
	    if(relres < tol)
	      y = H(1:i,1:i)\s(1:i);
          x = x+V(:,1:i)*y;
	      break;
        end
     
      end

      if(relres < tol)
        resvec(iter) = relres;
        break
      end
      y = H(1:m,1:m)\s(1:m);
      x = x+V(:,1:m)*y;                            
    
      % Update approximation
    
      r = b-A(x,varargin{:});                              
    
      % Compute residual
    
      s(i+1) = norm(r);
      relres = s(i+1)/bnrm;
    
      % Check convergence
    
      resvec(iter) = relres;
      if(relres < tol)
        break
      end
    
    end

  else
      
    % Stiffness matrix is given in terms of a double matrix  
    
    resvec = zeros(1,maxit);
    for iter = 1:maxit
      
      r = b-A*x;
      V(:,1) = r/norm(r);
      s = norm(r)*e1;
    
      % Construct orthonormal basis using Gram-Schmidt
    
      for i = 1:m
        w = A*V(:,i);
	    for k = 1:i,
	      H(k,i)= w'*V(:,k);
	      w = w-H(k,i)*V(:,k);
        end
	    H(i+1,i) = norm(w);
	    V(:,i+1) = w/H(i+1,i);
     
        % Apply Givens rotations
     
	    for k = 1:i-1
          temp = cs(k)*H(k,i)+sn(k)*H(k+1,i);
          H(k+1,i) = -sn(k)*H(k,i)+cs(k)*H(k+1,i);
          H(k,i) = temp;
        end
	    [cs(i),sn(i)] = rotmat(H(i,i),H(i+1,i)); 
      
        % Form i-th rotation matrix
      
        temp = cs(i)*s(i);                       
      
        % Approximate residual norm
      
        s(i+1) = -sn(i)*s(i);
	    s(i) = temp;
        H(i,i) = cs(i)*H(i,i)+sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
	    relres = abs(s(i+1))/bnrm;
     
        % Udpate approximation and exit
     
	    if(relres < tol)
	      y = H(1:i,1:i)\s(1:i);
          x = x+V(:,1:i)*y;
	      break;
        end
     
      end

      if(relres < tol)
        resvec(iter) = relres;
        break
      end
      y = H(1:m,1:m)\s(1:m);
      x = x+V(:,1:m)*y;                            
    
      % Update approximation
    
      r = b-A*x;                              
    
      % Compute residual
    
      s(i+1) = norm(r);
      relres = s(i+1)/bnrm;
    
      % Check convergence
    
      resvec(iter) = relres;
      if(relres < tol)
        break
      end
    
    end

  end
  
  % Check convergence
  
  if(relres > tol)
    flag = 0;
  else
    flag = 1;  
  end

  % Assign output arguments
  
  if(nargout == 1)
    varargout{1} = x;
  elseif(nargout == 2)
    varargout{1} = x
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
    varargout{5} = resvec(1:maxit);
  end
  
return

function [c,s] = rotmat(a,b)
% ROTMAT Givens rotation matrix parameters.
%
%   [C,S] = ROTMAT(A,B) computes the entries of the Givens rotation matrix
%   for A and B.
%
%   Example:
%
%   [c,s] = rotmat(a,b);

%   Iterative template routine --
%   Univ. of Tennessee and Oak Ridge National Laboratory October 1, 1993
%   Details of this algorithm are described in
%   "Templates for the Solution of Linear Systems: Building Blocks for
%   Iterative Methods",
%   Barrett, Berry, Chan, Demmel, Donato, Dongarra, Eijkhout, Pozo, Romine,
%   and van der Vorst,
%   SIAM Publications, 1993.
%   (ftp netlib2.cs.utk.edu; cd linalg; get templates.ps).

  if(b == 0.0)
    c = 1.0;
    s = 0.0;
  elseif(abs(b) > abs(a))
    temp = a/b;
    s = 1.0/sqrt(1.0+temp^2);
    c = temp*s;
  else
    temp = b/a;
    c = 1.0/sqrt(1.0+temp^2);
    s = temp*c;
  end

return