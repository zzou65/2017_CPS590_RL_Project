function [x,conv] = nmg_solve(mg_data,tol,maxit)
%NMG_SOLVE nested geometric multigrid solver
%
%   X = NMG_SOLVE(MG_DATA,TOL,MAXIT) computes the finite element solution
%   X of the Dirichlet problem described by the input argument MG_DATA
%   using the nested (geometric) multigrid method.
%
%   [X,CONV] = NMG_SOLVE(...) also returns a data structure CONV containing
%   information on the convergence of the iteration; see below.
%
%   The tolerance TOL and the maximal number of iterations MAXIT are used
%   as stopping criteria for the multigrid iterations on the various
%   levels.  These arguments may be either scalars or vectors of length
%   length(MG_DATA)-1.
%
%   If TOL is a vector, then TOL(i) is used as the tolerance on level i+1.
%   (Note that no multigrid iteration is performed on the coarsest level 
%   1.)  If TOL is a scalar, then the tolerance on level j is TOL*4^(LVL-j)
%   if LVL is the total number of levels.  In particular, the tolerance on
%   the finest level is TOL.
%
%   If MAXIT is a vector, then at most MAXIT(i) multigrid iterations are
%   done on level i+1.  If MAXIT is a scalar, then this value is used on
%   all levels.
%
%   Note that the solution is automatically calculated to a much higher
%   precision if MG_DATA is set to use relative errors.  This should be
%   avoided.
%
%   The output argument CONV is a struct containing informatiäon on the
%   convergence of the nested multigrid iteration.  Its fields are
%   identical to those of the analogous return argument of MG; they
%   essentially collect the data from the latter.  The fields are, if LVL
%   is the total number of levels, ie. LVL=length(MG_DATA) :
%     flag : a 1-by-LVL-1 logical array; CONV.flag(i) is true if the
%         multigrid iteration on level i+1 converged to the required
%         tolerance TOL(i) within MAXIT(i) steps.
%     iter : a 1-by-LVL-1 array; CONV.iter(i) is the number of computed
%         iterations on level i+1.
%     error : a 1-by-LVL-1 cell array containing estimated errors of the
%         multigrid iterations.
%     tol : a 1-by-LVL-1 array containing the tolerance at each level.
%     maxit : a 1-by-LVL-1 array containing the maximum number of
%         iterations on each level.
%     time.smooth : a 1-by-LVL-1 array containing cpu times required 
%         by the smoothing interations.
%     time.transfer : a 1-by-LVL-1 array containing cpu times required 
%         by the intergrid transfer operations.
%     time.coarse_solve : a 1-by-LVL-1 array containing cpu times required 
%         by the exact solver on the coarsest level during multigrid
%         iterations.
%     time.coarse_solve0 : time required by initial exact solver on
%         coarsest level.
%
%   For details on the data structure CONV, see the multigrid solver MG.
%   For details on the multigrid data structure MG_DATA, see its
%   initialization routines MG_MESH, MG_STIMA, MG_SMOOTH and MG_ERROR.
%   MG_DATA must be processed by all of these for the code to run.
%
%   Example :
%
%   Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%   mg_data = mg_mesh('mesh',Mesh,'ref',[2 7]);
%   mg_data = mg_stima(mg_data,...
%                     'f',@(x,varargin) -4*ones(size(x,1),1),...
%                     'gd',@(x,varargin) x(:,1).^2+x(:,2).^2);
%   mg_data = mg_smooth(mg_data,'m',[1 1],'sym_per',true);
%   mg_data = mg_error(mg_data,'rel',false);
%   u = mg_data{end}.u_bd;
%   [u(mg_data{end}.dofs),conv] = nmg_solve(mg_data,1e-6,20);
%
%   See also
%   initialization routines :
%     mg_mesh, mg_stima, mg_smooth, mg_error
%   other multigrid implementations :
%     mg, amg_solve, locmg_solve

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % initialize constants
  LVL = length(mg_data);
  
  if(numel(tol)==1)
    tol = tol*mg_data{end}.error_ctrl_scale.^(LVL-2:-1:0);
  end
  if(numel(maxit)==1)
    maxit = maxit*ones(1,LVL-1);
  end
  
  conv.flag = false(1,LVL-1);
  conv.iter = nan(1,LVL-1);
  conv.error = cell(1,LVL-1);
  conv.tol = tol;
  conv.maxit = maxit;
  conv.time.smooth = nan(1,LVL-1);
  conv.time.transfer = nan(1,LVL-1);
  conv.time.coarse_solve = nan(1,LVL-1);
  
  % initial guess on coarsest mesh
  t = cputime;
  if(mg_data{1}.n.free==0)
      x = zeros(0,1);
  else
      x = mg_data{1}.A\mg_data{1}.b;
  end
  conv.time.coarse_solve0 = cputime - t;
  
  % corrections on finer meshes
  for i = 1:LVL-1
    
    % prolong coarse grid solution to fine grid
    t = cputime;
    u = mg_data{i}.u_bd;
    u(mg_data{i}.dofs) = x;
    x = mg_data{i+1}.P_full(mg_data{i+1}.dofs,:)*u;
    time_prolong = cputime - t;
    
    % carry out multigrid iteration on level i+1
    [x,conv_step] = mg(x,{mg_data{1:i+1}},tol(i),maxit(i));
    
    % copy convergence information
    conv.flag(i) = conv_step.flag;
    conv.iter(i) = conv_step.iter;
    conv.error{i} = conv_step.error;
    conv.time.smooth(i) = conv_step.time.smooth;
    conv.time.transfer(i) = conv_step.time.transfer + time_prolong;
    conv.time.coarse_solve(i) = conv_step.time.coarse_solve;

  end

return