function varargout = mg_smooth(mg_data,varargin)
%MG_SMOOTH initialize smoothers for multigrid
%   
%   MG_DATA = MG_SMOOTH(MG_DATA,OPT) adds smoothing information to the
%   multigrid data structure MG_DATA.
%
%   OPT = MG_SMOOTH() returns the default arguments wrapped in a structure.
%
%   [MG_DATA,OPT] = MG_SMOOTH(MG_DATA,OPT) returns both structures.
%
%   The options can be given in one of two forms: either as a data
%   structure similar to the return argument OPT, or by specifying
%   parameters and values, see the examples below.  In the latter
%   case, if MG_DATA already contains smoother information, it is possible
%   to change only individual options, such as the number of smoothing
%   steps or the number of coarse grid cycles.
%
%   The option parameters are
%
%     m : number of smoothing steps, [default : [1 0] ].
%       m(1) : number of presmoothing steps.
%       m(end) : number of postsmoothing steps.
%
%     cyc : number of coarse grid corrections, [default : 1].  A value of 1
%         corresponds to V-cycles, 2 corresponds to W-cycles.
%
%     smoother : function handle for presmoother, [default : @gs_smooth].
%         The smoother must take at least the arguments x, A, b in that
%         order to carry out one iteration for A\b with initial guess x.
%
%     postsmoother : specify postsmoother, values :
%       'pre' [default] : use smoother as postsmoother.
%       (function handle) : use different smoother, see smoother.
%
%     args : additional arguments for smoother, values :
%       (cell array) : pass contents of args to smoother.
%       (function handle) : pass output of function handle applied to the
%         stiffness matrix to the smoother.
%
%     postargs : additional arguments for postsmoother, if postsmoother is
%         not 'pre'; see args.
%
%     per : permutation of indices for smoothing
%         To use different values for different smoothing steps, use a cell
%         array containing the desired values; these can be :
%       'CF' [default] : smooth coarse vertices, then fine vertices.
%       'FC' : smooth fine vertices, then coarse vertices.
%       'rand' : use random permutation.
%       'sort' : sort in ascending order according to per_fn.
%       'fn' : use permutation returned by per_fn.
%
%     postper : permutation of indices in postsmoothing, values :
%       'pre' [default] : use same permutation as for presmoothing.
%       (valid value of per argument) : see per.
%
%     per_fn : function handle used to determine permutation if per is
%         'sort' or 'fn'.  The function handle must take as arguments the
%         mesh and a logical array indexing the smoothed degrees of
%         freedom.  If per is 'sort', it should return an array of length
%         equal to the number of smoothed degrees of freedom; these are
%         smoothed in ascending order of the values in the array.  If per
%         is 'fn', the function should return an array of length equal to
%         the number of smoothed degrees of freedom containing the indices
%         in the order in which they should be smoothed.  To specify
%         different smoothing orders for different steps, wrap the function
%         handles in a cell array; then the function handles are iterated
%         in a circle.
%
%     per_postfn : same as per_fn for postsmoother.
%
%     sym_per : use opposite order for postsmoothing, [default : false].
%         If this argument is true and per is 'CF' or 'FC' and m is
%         symmetric and postper is 'pre', then the multigrid cycle with
%         Gauss-Seidel smoother is symmetric.
%
%     type : multigrid type, values :
%       'full' : smooth all vertices on each level.
%       'hier' : smooth only fine vertices (ie. hierarchical multigrid).
%
%   The multigrid data structure MG_DATA is a 1-by-LVL cell array of
%   structures, where LVL is the number of levels in multigrid.  For levels
%   lvl>1, MG_DATA{lvl} contains at least the following fields :
%
%     cyc : number of coarse grid corrections.
%
%     n.smooth : number of smoothed degrees of freedom (ie. number of
%         degrees of freedom included in smoothing times the number of
%         smoothing steps).
%
%     pre : presmoothing information.
%
%     post : postsmoothing information.
%
%   The structures MG_DATA{lvl}.pre and MG_DATA{lvl}.post contain at least
%   the following fields :
%
%     m : number of smoothing steps.
%
%     smoother : function handle for smoother.
%
%     dofs : degrees of freedom to smooth.
%
%     all : boolean specifying whether or not all degrees of freedom are
%         smoothed.
%
%     args : cell array containg arguments for smoother.
%
%     per : 1-by-m cell array containg permutations of all vertices to
%         apply before smoothing.
%
%     per_is_id : 1-by-m logical array specifying whether or not the
%         corresponding permutation is the identity.
%
%   Also, the options are stored in MG_DATA{1}.opt.smooth.
%
%   Examples :
%
%   1)  edit options directly
%
%       opt = mg_smooth;
%       opt.m = 1;
%       opt.sym_per = true;
%       mg_data = mg_smooth(mg_data,opt);
%
%   2)  or set individual parameters
%
%       mg_data = mg_smooth(mg_data,'m',1,'sym_per',true);
%
%   3)  change individual parameters (after (1) or (2))
%
%       mg_data = mg_smooth(mg_data,'m',2);
%
%   See also mg_mesh, mg_stima, mg_error, mg.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % return default options if no arguments are given
  if(nargin==0)
    smooth_opt.m = [1 0];
    smooth_opt.cyc = 1;
    
    smooth_opt.smoother = @gs_smooth;
    smooth_opt.postsmoother = 'pre'; % or some function handle like smoother
    smooth_opt.args = {}; % or cell or function handle with argument A
    smooth_opt.postargs = {}; % dito, overridden if postsmoother=='pre'
    
    smooth_opt.per = 'CF'; % permutation: 'CF','FC','rand','sort','fn'
    smooth_opt.postper = 'pre'; % permutation for postsmoothing
    smooth_opt.per_fn = []; % function handler for per
    smooth_opt.per_postfn = []; % function handler for postsmoothing per
    smooth_opt.sym_per = false; % use opposite order for postsmoothing
    
    smooth_opt.type = 'full'; % or 'hier' for hierarchical multigrid
    
    varargout{1} = smooth_opt;
    
    return
  end
  
  % read options
  if(length(varargin)==1)
    smooth_opt = varargin{1};
  else
    if(isfield(mg_data{1},'opt') && isfield(mg_data{1}.opt,'smooth'))
      smooth_opt = mg_data{1}.opt.smooth;
    else
      smooth_opt = mg_smooth;
    end
    for i=1:floor(0.5*length(varargin))
      smooth_opt.(varargin{2*i-1}) = varargin{2*i};
    end
  end
  
  % return options as second argument
  if(nargout>1)
    varargout{2} = smooth_opt;
  end
  
  % add options to first level of mg_data
  mg_data{1}.opt.smooth = smooth_opt;
  
  % preprocess smoothing options
  
  % postsmoother is same as presmoother (overrides options for
  % postsmoother)
  post_pre = isequal(smooth_opt.postsmoother,'pre');
  
  % type of permutation for smoothing
  if(isa(smooth_opt.per,'cell'))
    pre_per = smooth_opt.per;
  else
    pre_per = repmat({smooth_opt.per},1,smooth_opt.m(1));
  end
  
  if(isequal(smooth_opt.postper,'pre'))
    post_per0 = smooth_opt.per;
  else
    post_per0 = smooth_opt.postper;
  end
  
  if(isa(post_per0,'cell'))
    post_per = post_per0;
  else
    post_per = repmat({post_per0},1,smooth_opt.m(end));
  end
  
  % if type of permutation is 'sort' or 'fn', the followung function
  % handles are used
  pre_per_fn = smooth_opt.per_fn;
  post_per_fn = smooth_opt.per_postfn;
  if(isa(pre_per_fn,'function_handle'))
    pre_per_fn = {pre_per_fn};
  end
  if(post_pre)
    post_per_fn = pre_per_fn;
  elseif(isa(post_per_fn,'function_handle'))
    post_per_fn = {post_per_fn};    
  end
  n_pre_per_fn = numel(pre_per_fn);
  n_post_per_fn = numel(post_per_fn);
    
  
  % define smoothers on all levels
  for lvl=2:length(mg_data)
    
    % elementary parameters
    mg_data{lvl}.pre.m = smooth_opt.m(1);
    mg_data{lvl}.post.m = smooth_opt.m(end);
    mg_data{lvl}.cyc = smooth_opt.cyc;
    
    % function handles for smoothers
    mg_data{lvl}.pre.smoother = smooth_opt.smoother;
    if(isequal(smooth_opt.postsmoother,'pre'))
      mg_data{lvl}.post.smoother = smooth_opt.smoother;
    else
      mg_data{lvl}.post.smoother = smooth_opt.postsmoother;
    end
    
    % vertices to smooth
    if(isequal(smooth_opt.type,'full'))
      mg_data{lvl}.pre.dofs = true(mg_data{lvl}.n.free,1);
      mg_data{lvl}.post.dofs = mg_data{lvl}.pre.dofs;
    elseif(isequal(smooth_opt.type,'hier'))
      mg_data{lvl}.pre.dofs = true(mg_data{lvl}.n.free,1);
      mg_data{lvl}.pre.dofs(1:mg_data{lvl-1}.n.free) = false;
      mg_data{lvl}.post.dofs = mg_data{lvl}.pre.dofs;
    else
      error('Unknown smoothing pattern');
    end
    mg_data{lvl}.pre.all = all(mg_data{lvl}.pre.dofs);
    mg_data{lvl}.post.all = all(mg_data{lvl}.post.dofs);
    n_smooth = nnz(mg_data{lvl}.pre.dofs);
    mg_data{lvl}.n.smooth = n_smooth*(mg_data{lvl}.pre.m+mg_data{lvl}.post.m);
    
    % presmoother arguments
    if(isa(smooth_opt.args,'function_handle'))
      mg_data{lvl}.pre.args = {smooth_opt.args(...
        mg_data{lvl}.A(mg_data{lvl}.pre.dofs,mg_data{lvl}.pre.dofs))};
    else
      mg_data{lvl}.pre.args = smooth_opt.args;
    end
    
    % postsmoother arguments
    if(post_pre)
      mg_data{lvl}.post.args = mg_data{lvl}.pre.args;
    else
      if(isa(smooth_opt.postargs,'function_handle'))
        mg_data{lvl}.post.args = {smooth_opt.postargs(...
          mg_data{lvl}.A(mg_data{lvl}.post.dofs,mg_data{lvl}.post.dofs))};
      else
        mg_data{lvl}.post.args = smooth_opt.postargs;
      end
    end
    
    % permutation of indices in presmoother
    per = find(mg_data{lvl}.pre.dofs);
    dofs = mg_data{lvl}.dofs;
    dofs(dofs) = mg_data{lvl}.pre.dofs;
    mg_data{lvl}.pre.per_is_id = false(1,mg_data{lvl}.pre.m);
    for i=1:mg_data{lvl}.pre.m
      switch pre_per{i}
        case 'CF'
          mg_data{lvl}.pre.per{i} = per;
        case 'FC'
          mg_data{lvl}.pre.per{i} = per(end:-1:1);
        case 'rand'
          mg_data{lvl}.pre.per{i} = per(randperm(length(per)));
        case 'sort'
          [dummy,ind_per] = ...
            sort(pre_per_fn{mod(i-1,n_pre_per_fn)+1}(mg_data{lvl}.mesh,dofs));
          mg_data{lvl}.pre.per{i} = per(ind_per);
        case 'fn'
          ind_per = pre_per_fn{mod(i-1,n_pre_per_fn)+1}(mg_data{lvl}.mesh,dofs);
          mg_data{lvl}.pre.per{i} = per(ind_per);
      end
      mg_data{lvl}.pre.per_is_id(i) = isequal(per,mg_data{lvl}.pre.per{i});
    end
          
    % permutation of indices in postsmoother
    per = find(mg_data{lvl}.post.dofs);
    dofs = mg_data{lvl}.dofs;
    dofs(dofs) = mg_data{lvl}.post.dofs;
    mg_data{lvl}.post.per_is_id = false(1,mg_data{lvl}.post.m);
    for i=1:mg_data{lvl}.post.m
      switch post_per{i}
        case 'CF'
          mg_data{lvl}.post.per{i} = per;
        case 'FC'
          mg_data{lvl}.post.per{i} = per(end:-1:1);
        case 'rand'
          mg_data{lvl}.post.per{i} = per(randperm(length(per)));
        case 'sort'
          [dummy,ind_per] = ...
            sort(post_per_fn{mod(i-1,n_post_per_fn)+1}(mg_data{lvl}.mesh,dofs));
          mg_data{lvl}.post.per{i} = per(ind_per);
        case 'fn'
          ind_per = post_per_fn{mod(i-1,n_post_per_fn)+1}(mg_data{lvl}.mesh,dofs);
          mg_data{lvl}.post.per{i} = per(ind_per);
      end
      if(smooth_opt.sym_per)
        mg_data{lvl}.post.per{i} = mg_data{lvl}.post.per{i}(end:-1:1);
      end
      mg_data{lvl}.post.per_is_id(i) = isequal(per,mg_data{lvl}.post.per{i});
    end
    
  end % for loop over levels, index lvl
  
  % return mg_data multigrid data structure as first argument
  varargout{1} = mg_data;
  
return