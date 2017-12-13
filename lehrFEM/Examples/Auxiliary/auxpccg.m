function [itnum,flag,relres,residuals] = auxpccg(Mesh,ALPHA,BETA,tol,maxit)
% Solve H(curl)-elliptic problem by means of a CG iteration
% Use constant vectorfield (1,1) as source function and 0 as initial guess
%
% Nodal auxiliary space preconditioning with zero Dirichlet boundary conditions
% imposed on auxiliary space. 
%
% ALPHA -> handles to functions of type @(x,varargin), where
% BETA     the rows of x contain positions 
% (ALPHA,BETA can also be mere positive numbers)
%
% tol -> termination of PCG: relative decrease of residual
% maxit -> maximal number of CG steps. 
%
% RETURN VALUES: (see doucmentation of pcg_solve)
% itnum -> number of PCG steps
% flag  -> the integer FLAG
%   specifies wheter the method has converged or not:
%    1 Method has converged within MAXIT steps to the prescribed tolerance
%      TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%

if (nargin < 5), maxit = 200; end
if (nargin < 4), tol = 1.0E-6; end
if (nargin < 3), BETA = 1; end
if (nargin < 2), ALPHA = 1; end
if (nargin == 0), error('mesh argument required'); end

% Constant source vectorfield

F = @(x,varargin) ones(size(x,1),2);

% Obtain matrices 

[A,P,VL,Grad,Lpot,f] = auxpcops(Mesh,ALPHA,BETA,F);

% fprintf('size(A) = '); disp(size(A));
% fprintf('size(f) = '); disp(size(f));

% Preconditioner defined through local function below
% Select Jacobi smoother for preconditioner (last argument = 0 )
% Select Gauss-Seidel smoother (last argument > 0)
B = @(r) curlauxpc(r,A,P,VL,Grad,Lpot,2);

% Invoke PCG iteration 
x = zeros(size(A,1),1);
[x,flag,relres,itnum,residuals] = pcg_solve(x,A,f,tol,maxit,B);



