function [x,flag,relres,iter,resvec,times] = namg_solve(b,A,AMGOptions,tol,maxit)
%NAMG_SOLVE nested algebraic multigrid solver
%
%   A = NAMG_SOLVE(B,A,AMGOPTIONS,TOL,MAXIT)
%   computes the solution X = A\B using the nested algebraic multigrid
%   method. The tolerance TOL and the maximal number of iterations MAXIT
%   are used as stopping criteria for the multigrid iterations.
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
%   See AMG_SOLVE for a list of parameters.
%
%   [X,FLAG] = NAMG_SOLVE(...) the integer FLAG specifies whether or
%   not the method has converged:
%    1 Method converged within MAXIT steps to the prescribed tolerance TOL.
%    0 Method did not converge within MAXIT steps to the prescribed
%      tolerance TOL.
%
%   [X,FLAG,RELRES] = NAMG_SOLVE(...) also returns the relative
%   residual RELRES at the last multigrid iteration, i.e. RELRES is the
%   norm of x_k - x_{k-1} diveded by the norm of the right hand side if the
%   last iteration is k.
% 
%   [X,FLAG,RELRES,ITER] = NAMG_SOLVE(...) also returns the number of 
%   multigrid iterations ITER that were performed on each level.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = NAMG_SOLVE(...) also returns the relatice
%   residual RESVEC at every iteration on the highest level on the nested
%   iteration.
%
%   [X,FLAG,RELRES,ITER,RESVEC,TIMES] = NAMG_SOLVE(...) also returns a
%   structure containing the times required for various parts of the code. 
%   In particular, 'TIMES.level{lvl}.V_cycles' contains information on the
%   time dedicated to V-cycles on level lvl of the nested iteration.
%   
%   See also amg_solve, nmg_solve.

%   Copyright 2007-2007 Claude Gittelson
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
  end
  
  % Initialize constants
  
  LVL = length(AMGData);
  
  tol_C = tol*4^(LVL-2);            % tolerance on second coarsest grid
  
  flag = 1;
  relres = zeros(1,LVL-1);
  iter = zeros(1,LVL-1);
  
  % Calculate right hand side on coarse grids
  
  B = cell(1,LVL);
  B{1} = b;
  for lvl = 2:LVL
    B{lvl} = AMGData{lvl-1}.R*B{lvl-1};
  end
  
  % Initial guess on coarsest mesh

  if(isempty(AMGData{LVL}.A))
      x = zeros(0,1);
  else
      x = AMGData{LVL}.A\B{LVL};
  end
  
  % corrections on finer meshes

  for lvl = 1:LVL-1

    % carry out multigrid iterations on level lvl

    x = AMGData{LVL-lvl}.P*x;
    [x,flag1,relres(lvl),iter(lvl),resvec,times.level{lvl}] = ...
      amg_solve(x,B{LVL-lvl},[],{AMGData{(LVL-lvl):LVL}},tol_C,maxit);
    flag = flag && flag1;

    % update error tolerance

    tol_C = tol_C/4;

  end

return
  
  