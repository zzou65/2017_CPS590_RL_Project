function [x,flag,relres,iter,resvec,times,AMGData] = amg_solve(x,b,A,AMGOptions,tol,maxit)
%AMG_SOLVE algebraic multigrid solver (wrapper)
%   
%   X = AMG_SOLVE(X0,B,A,AMGOPTIONS,TOL,MAXIT)
%   computes the solution X = A\B from the initial guess X0 using the
%   algebraic multigrid method.  The tolerance TOL and the maximal number
%   of iterations MAXIT are used as stopping criteria for the multigrid
%   iterations.
%
%   AMGOPTIONS is a structure containing multigrid options.  If AMGOptions
%   is empty, the default value is used,
%       AMGOPTIONS = AMGDefaultOptions;
%   To use different options, load the default options,
%       myAMGOptions = AMGDefaultOptions;
%   and change the desired values.  See below for a list of parameters.
%   Alternatively, this argument can be replaced by a cell array AMGDATA
%   containing multilevel information, see below.  AMGDATA is constructed
%   from AMGOPTIONS by
%       AMGDATA = AMGSetup(A,AMGOPTIONS);
%   where A is the stiffness matrix.  AMGOPTIONS can be recovered from
%   AMGDATA by
%       AMGOPTIONS = AMGDATA{1}.opt;
%   If the AMGDATA argument is given, A is ignored.
%
%   [X,FLAG] = AMG_SOLVE(...) the integer FLAG specifies whether or
%   not the method has converged:
%    1 Method converged within MAXIT steps to the prescribed tolerance TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%
%   [X,FLAG,RELRES] = AMG_SOLVE(...) also returns the relative
%   residual RELRES at the last multigrid iteration, i.e. RELRES is the
%   norm of x_k - x_{k-1} diveded by the norm of the right hand side if the
%   last iteration is k.
% 
%   [X,FLAG,RELRES,ITER] = AMG_SOLVE(...) also returns the number of 
%   multigrid iterations ITER that were performed.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = AMG_SOLVE(...) also returns the
%   values of the relative residuals RESVEC at each multigrid iteration.
%
%   [X,FLAG,RELRES,ITER,RESVEC,TIMES] = AMG_SOLVE(...) also returns a
%   structure containg the times required for various parts of the code.
%   In particular, 'TIMES.V_cycles' contains onfromation on the time
%   dedicated to V-cycles.
%
%   [X,FLAG,RELRES,ITER,RESVEC,TIMES,AMGDATA] = AMG_SOLVE(...) also returns
%   a cell array AMGDATA of structures containing multigrid information on
%   all the levels.  This can be reused as an input argument in place of
%   AMGOPTIONS to significantly speed up the code.
%
%
%   LIST OF PARAMETERS
%   These are some of the fields of AMGOPTIONS.
%
%   Note : this code always uses V-cycles.
%
%
%   LEVELS
%
%   .levsmax : maximal number of levels [default : 0, ie. no limit]
%       Setting .levsmax to 1 generates an error.
%
%   .mincoarse : minimal number of dofs for coarsening [default : 100]
%       If there are fewer grid points, then no more coarsening is done and
%       an exact solver is used on the current level.
%
%   .max_coarse_frac : maximal relative number of coarse grid points
%       [default : 0.75]
%       Setting this option to small often generates errors.
%
%
%   SMOOTHING
%
%   .pre.its : number of presmoothing iterations [default : 2]
%
%   .post.its : number of postsmoothing iterations [default : 2]
%
%   .pre.type : type of presmoother; values :
%     'Gauss-Seidel CF' [default] : Gauss-Seidel smoother
%         Additional parameters : 
%         .pre.GSperm : ordering of indices in Gauss-Seidel
%             [default : [], no special ordering]
%             argument should be an array containing the indices in the
%             order in which they should be smoothed.
%     'Gauss-Seidel upper CF' : "transposed" Gauss-Seidel smoother
%     'Jacobi' : Jacobi smoother
%         Additional parameters :
%         .pre.jac_omega : smoother is pre.jac_omega*(Jacobi smoother)
%             [default : 0, ie use upper bound for eigenvalues of D^{-1}*A]
%     'ainv' : approximate factorization of inverse of SPD matrix
%         Ainv = Z*inv(D)*Z.' for Z upper triangular, D diagonal
%         Additional parameters :
%         .pre.ainv_tol : drop tolerance in ainv routine [default : 0.1]
%         .pre.ainv_omega : smoother is .pre.ainv_omega*Ainv [default : 1]
%         .pre.ainv_shift : ainv is applied to A + .pre.ainv_shift*I
%             [default : 0.05]
%
%   .post.pre : use presmoother as postsmoother [default : 1]
%
%   .post.type : if .post.pre==0, type of postsmoother; same values as .pre.type.
%
%
%   COARSE GRID SELECTION
%
%   .CF.theta : sensitivity of strong connectedness [default : 0.25]
%       Vertex i is strongly connected to j if 
%         -A(i,j) >= .CF.theta*max(-A(i,:))                            (sc)
%       where diagonal entries are not included in the maximum on the right
%       hand side.  In the following, As is the matrix containing the
%       off-diagonal strong entries of A, ie A(i,j) is in As(i,j) iff j is
%       strongly connected to i.
%
%   .CF.use_abs : use absolute value of A instead of -A in (sc)
%       [default : 0]
%
%   .CF.use_diag : use diagonal entries in (sc) and in As [default : 0]
%
%   .CF.direction : direction in (sc); values : 
%     'row' [default] : as in (sc)
%     'row and column' : i is strongly connected to j if (sc) holds for A
%         or for the transposed of A.
%
%   .CF.make_unconnected_F : make vertices with no strong connections (in
%       both directions) fine grid points [default : 1]
%
%   .CF.check_F_F : find strongly connected fine grid points with no common
%       coarse grid point; values :
%     0 : do nothing
%     1 [default] : check rows of As and add necessary coarse grid points
%     2 : check rows and columns of As and add necessary coarse grid points
%
%
%   INTERPOLATION
%
%   .P.method : interpolation method
%       The fine-grid entries of the prolongation matrix approximate
%         -A(f,f)\A(f,c)
%       where f denotes the fine grid indices and c the coarse indices.
%       The restriction matrix is the transpose of the prolongation matrix.
%       Values :
%     'Ruge-Stueben' [default], 'R-S' : Ruge-Stueben interpolation
%         Additional parameters : see below
%     'ainv' : approximate factorization of inverse of SPD matrix
%         Ainv = Z*inv(D)*Z.' for Z upper triangular, D diagonal
%         Additional parameters :
%         .P.ainv_tol : drop tolerance in ainv routine [default : 0.1]
%     'ainv strong' : ainv for As instead of A (only strong connections)
%     'gsai' : generalised sparse approximate inverse
%
%   .P.droptol : drop tolerance for prolongation matrix [default : 0.0]
%
%
%   RUGE-STUEBEN INTERPOLATION
%
%   .P.method = 'Ruge-Stueben' [default] or (equivalent) .P.method = 'R-S'
%
%   .P.method64 : variations of Ruge-Stueben interpolation formula; only
%       strong connections are used; values :
%     'original' [default] : original method
%     'original abs' : abs in sum
%     'original abs abs' : abs in sum and numerator
%     'new1' : different interpolaiton formula
%     'new2' : variation of different interpolation formula
%
%   .P.weak : handling of weak connections; values :
%     'lump' [default] : lump values in the diagonal
%     'scale' : ignore weak connections and rescale weights for strong
%         connections
%     'trash' : ignore weak connections
%
%   .P.nofill64 : prevent fillup in Ruge-Stueben interpolation
%       [default : 1]
%       Only strong connections are used in interpolation.
%
%
%   Examples:
%
%   1)  default options
%       
%       x = amg_solve(zeros(size(b)),b,A,[],1e-6,20);
%
%   2)  only one presmoothing step and one postsmoothing step
%
%       myAMGOptions = AMGDefaultOptions;
%       myAMGOptions.pre.its = 1;
%       myAMGOptions.post.its = 1;
%       x = amg_solve(zeros(size(b)),b,A,myAMGOptions,1e-6,20);
%
%   3)  define more connections as strong connections in coarse grid
%       selection
%
%       myAMGOptions = AMGDefaultOptions;
%       myAMGOptions.P.theta = 0.2;
%       x = amg_solve(zeros(size(b)),b,A,myAMGOptions,1e-6,20);
%
%
%   See also multigrid_solve, AMGSetup, AMGSelectCoarseGrid,
%   AMGMakeInterpolation, AMGDefaultOptions, Ruge_Stueben, ainv,
%   ainv_explicit, AMGVcycle, namg_solve

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Get input arguments
  
  if(isempty(AMGOptions))
    AMGOptions = AMGDefaultOptions;
  end
  
  % Generate AMG data structure
  
  if(isa(AMGOptions,'struct'))
    [AMGData,times] = AMGSetup(A,AMGOptions);
  else
    AMGData = AMGOptions;
    A = AMGData{1}.A;
  end
  
  % Initialize timer
  
  times.V_cycles.elapsed = 0;
  times.V_cycles.nintervals = 0;

  % Initialize
  
  bnrm = norm(b);
  if(bnrm == 0)
    bnrm = 1;
  end
  y = b - A*x;
  
  % Run multigrid solver
  
  iter = 0;
  relres = tol+1;
  resvec = zeros(1,maxit);
  while(iter < maxit && relres > tol)
    t = cputime;
    c = AMGVcycle(AMGData,y);
    times.V_cycles.elapsed = times.V_cycles.elapsed + cputime - t;
    times.V_cycles.nintervals = times.V_cycles.nintervals + 1;
    
    y = y - A*c;
    x = x + c;
    
    iter = iter+1;
    relres = norm(c)/bnrm;
    resvec(iter) = relres;
  end

  % Check convergence
  
  if(relres > tol)
    flag = 0;
  else
    flag = 1;  
  end
  
  resvec = resvec(1:iter);
  
return