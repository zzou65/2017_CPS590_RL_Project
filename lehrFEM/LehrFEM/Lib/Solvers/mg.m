function [x,conv] = mg(varargin)
%MG geometric multigrid solver and preconditioner
%   
%   X = MG(...) computes the finite element solution X of the Dirichlet
%   problem described by the input arguments using the (unnested geometric)
%   multigrid method.
%
%   [X,CONV] = MG(...) also returns a data structure CONV containing
%   information on the convergence of the multigrid algorithm.
%
%   Valid syntaxes for the input arguments are
%     MG(X0,B,MG_DATA,TOL,MAXIT)
%     MG(X0,MG_DATA,TOL,MAXIT) if B is defined in MG_DATA
%     MG(MG_DATA,TOL,MAXIT) uses X0=0 as initial guess
%     MG(B,MAXIT,MG_DATA) uses X0=B as initial guess
%   Note that the last version can be used as a preconditioner for a Krylov
%   space method.
%
%   The solution vector X, the initial guess X0 and the load vector B are
%   all defined only on the non-Dirichlet-boundary vertices of the grid.
%   Inhomogeneous Dirichlet boundary conditions should be incorporated in
%   the load vector.  
%
%   The tolerance TOL and the maximum number of iterations MAXIT are used 
%   as stopping criteria for the multigrid iteration.  If TOL is not given,
%   or if no error norm is specified in MG_DATA (see MG_ERROR), then 
%   exactly MAXIT iterations are performed.  If an error norm is specified,
%   then the iteration stops when this value reaches the tolerance TOL.
%
%   The output argument CONV is a struct containing information on the
%   convergence of the multigrid iteration.  It has the following fields :
%     flag : true if the iteration converged to the required tolerance TOL
%         within MAXIT steps in the norm specified, false otherwise, or if
%         no norm is given.
%     iter : number of computed iterations.
%     error.(nrm) : a 1-by-iter vector containing the error in the norm nrm
%         after every multigrid cycle.
%     time.smooth : time required for smoothing.
%         .transfer : time required for intergrid transfer, ie. for
%             prolongation and restriction
%         .coarse_solve : time required for solving the equation exactly on
%             the coarsest grid.
%         .error : time required for calculating errors.
%
%   For details on the multigrid data structure MG_DATA, see its
%   initialization routines MG_MESH, MG_STIMA, MG_SMOOTH and MG_ERROR.
%   MG_DATA must be processed by all of these for the code to run as a
%   solver and by the first three for the code to run as a preconditioner
%   (ie. with the last syntax listed above).
%
%   Examples :
%
%   1)  MG as solver
%
%       Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%       mg_data = mg_mesh('mesh',Mesh,'ref',[2 7]);
%       mg_data = mg_stima(mg_data,...
%                         'f',@(x,varargin) -4*ones(size(x,1),1),...
%                         'gd',@(x,varargin) x(:,1).^2+x(:,2).^2);
%       mg_data = mg_smooth(mg_data,'m',[1 1],'sym_per',true);
%       mg_data = mg_error(mg_data,'energy',false,'l2',true,'h1',true);
%       u = mg_data{end}.u_bd;
%       [u(mg_data{end}.dofs),conv] = mg(mg_data,1e-6,20);
%
%   2)  MG as preconditioner
%
%       Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%       mg_data = mg_mesh('mesh',Mesh,'ref',[2 7]);
%       mg_data = mg_stima(mg_data,...
%               'f',@(x,varargin) -4*ones(size(x,1),1),...
%               'gd',@(x,varargin) x(:,1).^2+x(:,2).^2);
%       mg_data = mg_smooth(mg_data,'m',1,'sym_per',true);
% 
%       u = mg_data{end}.u_bd;
%       A = mg_data{end}.A;
%       b = mg_data{end}.b;
%       x0 = zeros(size(b));
%       u(mg_data{end}.dofs) = pcg_solve(x0,A,b,1e-6,20,@mg,1,mg_data);
%
%   See also
%   initialization routines :
%     mg_mesh, mg_stima, mg_smooth, mg_error
%   other multigrid implementations :
%     nmg_solve, amg_solve, locmg_solve

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize parameters
  error_flag = true;          % calculate / estimate errors
  error_ctrl_flag = true;     % control iteration by error estimate.  If false, maxit iterations are done
  conv_flag = (nargout>=2);   % store convergence data

  % read input arguments
  
  % find mg_data (only cell array)
  for mg_data_ind=1:nargin
    if(isa(varargin{mg_data_ind},'cell'))
      mg_data = varargin{mg_data_ind};
      LVL = length(mg_data);
      break;
    end
  end
  
  % read arguments and set parameters for preconditioner
  if(mg_data_ind==nargin && mg_data_ind==3)
    b = varargin{1};
    x = b;
    maxit = varargin{2};
    error_flag = false;
    
  % read arguments and set parameters for solver
  else
    tol = varargin{mg_data_ind+1};
    maxit = varargin{mg_data_ind+2};
    if(mg_data_ind>1)
      x = varargin{1};
    else
      x = zeros(mg_data{LVL}.n.free,1);
    end
    if(mg_data_ind>2)
      b = varargin{2};
    else
      b = mg_data{LVL}.b;
    end
    error_ctrl_flag = mg_data{LVL}.error_ctrl_flag;
    error_ctrl = mg_data{LVL}.error_ctrl;
  end
  
  % initialize convergence data structure conv
  if(conv_flag)
    
    % convergence of iteration
    conv.flag = false;
    
    % error estimators
    if(error_flag)
      error_names = fieldnames(mg_data{LVL}.error);
      for i=1:length(error_names)
        conv.error.(error_names{i}) = nan(1,maxit);
      end
    end
    
    % timer
    conv.time.smooth = 0;
    conv.time.transfer = 0;
    conv.time.coarse_solve = 0;
    conv.time.error = 0;
  end
  
  % run multigrid solver
  for iter=1:maxit
    
    % multigrid cycle
    x0 = x;
    [x,timer] = mg_step(x0,b,mg_data,LVL);
    
    % update global timers
    if(conv_flag)
      conv.time.smooth = conv.time.smooth + timer.smooth;
      conv.time.transfer = conv.time.transfer + timer.transfer;
      conv.time.coarse_solve = conv.time.coarse_solve + timer.coarse_solve;
    end
    
    % calculate error and check for convergence
    if(error_flag) % need to estimate error
      t = cputime;
      if(conv_flag) % ... in all norms and store values
        for i=1:length(error_names)
          conv.error.(error_names{i})(iter) = mg_data{LVL}.error.(error_names{i})(x,x0,mg_data{LVL}.A,b);
        end
        if(error_ctrl_flag)
          err = conv.error.(error_ctrl)(iter);
        end
        conv.time.error = conv.time.error + cputime - t;
      elseif(error_ctrl_flag) % ... in just one norm, discard value
        err = mg_data{LVL}.error.(error_ctrl)(x,x0,mg_data{LVL}.A,b);
      end
      if(error_ctrl_flag)
        if(mg_data{LVL}.error_rel && iter==1) % use relative errors
          tol = tol*err;
        end
        if(err<tol) % check for convergence
          if(conv_flag)
            conv.flag = true;
          end
          break;
        end
      end
    end
  end
  
  
  % define output arguments
  if(conv_flag)
    
    % number of iteraions
    conv.iter = iter;
    
    % errors in various norms
    if(error_flag)
      for i=1:length(error_names)
        conv.error.(error_names{i}) = conv.error.(error_names{i})(1:iter);
      end
      % check for convergence if iteration not controlled by error
      % estimator
      if(~error_ctrl_flag && ~isempty(error_names)) 
        conv.flag = true;
        for i=1:length(error_names)
          conv.flag = conv.flag & (conv.error.(error_names{i})(iter)<tol);
        end
      end
    end
  end
  
return
  
  
%%% multigrid cycle %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,timer] = mg_step(x,b,mg_data,lvl)
%MG_STEP geometric multigrid cycle
%
%   [X,TIMER] = MG_STEP(X0,B,MG_DATA,LVL) recursively calculates a multigrid
%   cycle with initial guess X0, load vector B and multigrid data MG_DATA
%   starting at level LVL.  The return arguments are the solution X and the
%   struct TIMER that stores the times required for smoothing, intergrid
%   transfer and direct solving on the coarsest mesh in the fields
%   'smooth', 'transfer' and 'coarse_solve', respectively.
%
%   See alse mg, mg_smooth.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize timer
  timer.smooth = 0;
  timer.transfer = 0;
  timer.coarse_solve = 0;

  % use direct solver on coarsest mesh
  if(lvl==1)
    t = cputime;
    if(mg_data{1}.n.free>0)
      x = mg_data{1}.U\(mg_data{1}.L\b);
    else
      x = zeros(0,1);
    end
    timer.coarse_solve = cputime - t;
    return
  end
  
  % initialize constants
  A = mg_data{lvl}.A;
  P = mg_data{lvl}.P;
  
  % presmooth
  [x,time_presmooth] = smooth(x,A,b,mg_data{lvl}.pre);
  
  % coarse grid correction
  if(mg_data{lvl}.cyc)
    t = cputime;
    res = P'*(b-A*x);
    time_restrict = cputime - t;
    cor = zeros(size(res));
    for i = 1:mg_data{lvl}.cyc
      [cor,ctimer] = mg_step(cor,res,mg_data,lvl-1);
      timer.smooth = timer.smooth + ctimer.smooth;
      timer.transfer = timer.transfer + ctimer.transfer;
      timer.coarse_solve = timer.coarse_solve + ctimer.coarse_solve;
    end
    t = cputime;
    x = x+P*cor;
    time_prolong = cputime - t;
  else
    time_restrict = 0;
    time_prolong = 0;
  end
  
  % postsmooth
  [x,time_postsmooth] = smooth(x,A,b,mg_data{lvl}.post);
  
  % update timer
  timer.smooth = timer.smooth + time_presmooth + time_postsmooth;
  timer.transfer = timer.transfer + time_restrict + time_prolong;
  
return


%%% smoothing iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,time_smooth] = smooth(x,A,b,data)
%SMOOTH apply multigrid smoother
%
%   [X,TIME] = SMOOTH(X0,A,B,DATA) performs DATA.M smoothing steps with
%   smoother DATA.smoother, initial guess X0, matrix A, right hand side B
%   and additional parameters specified by DATA.  The return arguments are
%   the solution X and the elapsed time TIME.
%
%   See also mg, mg_smooth.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  t = cputime;
  
  % determine right hand side
  if(~data.all)
    notdofs = ~data.dofs;
    b(data.dofs) = b(data.dofs) - A(data.dofs,notdofs)*x(notdofs);
  end
  
  % do smoothing
  if(all(data.per_is_id) && data.all) % use most efficient code
    for i = 1:data.m
      x = data.smoother(x,A,b,data.args{:});
    end
    
  elseif(all(data.per_is_id)) % use code with smoothing of only some of the nodes
    As = A(data.dofs,data.dofs);
    bs = b(data.dofs);
    y = x(data.dofs);
    for i = 1:data.m
      y = data.smoother(y,As,bs,data.args{:});
    end
    x(data.dofs) = y;
    
  else % use more flexible code with permutations of indices and smoothing of only some of the nodes
    for i = 1:data.m
      per = data.per{i};
      x(per) = data.smoother(x(per),A(per,per),b(per),data.args{:});
    end
  end
  
  time_smooth = cputime - t;
  
return