function [lmin,lmax,flag,iter] = auxpccond(Mesh,ALPHA,BETA,tol,maxit)
% Compute the condition number of the matrix related to bilinear form
% (\alpha\curl u,\curl v) + (\beta u,v) when preconditioned by means
% of the nodal auxiliary space preconditioner.
%
% * Zero Dirichlet boundary conditions everywhere
% * Zero boundary values imposed on auxiliary space
%
% Mesh -> complete mesh data structure with edges
%
% ALPHA -> handles to functions of type @(x,varargin), where
% BETA     the rows of x contain positions 
% (ALPHA,BETA can also be mere positive numbers)
%
% tol -> tolerance for Lanczos condition number estimate
% maxit -> maximal dimension of Krylov subspace for condition number
%          estimate
%
% RETURN VALUES: 
% 
% lmin,lmax : extremal eigenvalues of preconditioned system
% flag :   0 Lnczos converged within maxit steps to prescribed tolerance TOL.
%          1 Lanczos method failed to converge.
% iter : Krylov subspace dimension for Lanczos method

if (nargin < 5), maxit = 200; end
if (nargin < 4), tol = 1.0E-6; end
if (nargin < 3), BETA = 1; end
if (nargin < 2), ALPHA = 1; end

if (nargin == 0), error('mesh argument required'); end

% Obtain matrices 

[A,P,VL,Grad,Lpot] = auxpcops(Mesh,ALPHA,BETA);

% Preconditioner defined through local function below
% Select Jacobi smoother for preconditioner (last argument = 0)
% Select Gauss-Seidel smoother for preconditioner (last argument > 0)
B = @(r) curlauxpc(r,A,P,VL,Grad,Lpot,2);

[lmin,lmax,flag,iter] = eigen_Lanzcos(A,tol,maxit,B); 

end

% ######################################################################







