function varargout = mg_error(mg_data,varargin)
%MG_ERROR define error functionals for multigrid iteration
%   
%   MG_DATA = MG_ERROR(MG_DATA,OPT) adds information on error functionals
%   to monitor and control the multigrid iteration to the multigrid data
%   structure MG_DATA.
%
%   OPT = MG_ERROR() returns the default arguments wrapped in a structure.
%
%   [MG_DATA,OPT] = MG_ERROR(MG_DATA,OPT) returns both structures.
%
%   The options can be given in one of two forms: either as a data
%   structure similar to the return argument OPT, or by specifying
%   parameters and values, see the examples below.  In the latter
%   case, if MG_DATA already contains error information, it is possible
%   to change only individual options, such as the number of smoothing
%   steps or the number of coarse grid cycles.
%
%   The option parameters are
%
%     iter : estimate error of iterate x_n by norm of x_n - x_{n-1},
%         [default : true].
%
%     exact : calculate exact errors, [default : false].
%         If true, the exact solution is calculated and used to determine
%         the errors.  Note that, in particular, this requires that MG_DATA
%         contains information on the stiffness matrix and right hand side.
%
%     eucl : use euclidean norm, [default : false].
%
%     l2 : use L2 norm, [default : false].
%
%     h1 : use H1 seminorm, [default : false].
%
%     energy : use energy norm, [default : true].
%         Note that this requires information on the stiffness matrix in
%         MG_DATA.
%
%     custom : cell array of functon handles for custom error functionals.
%         The function handles should take the arguments x_n, x_{n-1}, A, b
%         in that order and return a scalar.
%
%     custom_names : names of custom error functionals.  If these are not
%         given, 'custom1', 'custom2', etc are used.
%
%     ctrl : name of error functional to use to control the behavior of the
%         multigrid algorithm, values :
%       'auto' [default] : use the first of the following list that is
%           being calculated (ie. to use 'l2_iter', l2 and iter must be
%           true).
%       'eucl_iter' : euclidean norm of x_n - x_{n-1}.
%       'l2_iter' : L2 norm of x_n - x_{n-1}.
%       'h1_iter' : H1 seminorm of x_n - x_{n-1}.
%       'energy_iter' : energy norm of x_n - x_{n-1}.
%       'eucl_exact' : euclidean norm of exact error vector.
%       'l2_exact' : L2 norm of exact error.
%       'h1_exact' : H1 seminorm of exact error.
%       'energy_exact' : energy norm of exact error.
%       (custom_names{i}) : i-th custom error functional.
%       '' : do not use any error estimator to control iteration.  Note
%           that error estimators may still be used to monitor the
%           convergence of a multigrid iteration.
%
%     ctrl_scale : scaling behaviour of discretization error wrt ctrl norm,
%         ie. a regular refinement decreases the discretization error by
%         this factor.  Should be 4 for L2-type norms and 2 for H1-type
%         norms. Values :
%       'auto' [default] use 4 for euclidean and L2 norms and 2 for H1 and
%           energy norms.
%       (scalar) : use this value
%
%     rel : use relative errors, [default : false].
%
%     add_scapro : add scalar products used in L2 and H1 norms to multigrid
%         data structure if these norms are used, [default : true].
%
%   The multigrid data structure MG_DATA is a 1-by-LVL cell array of
%   structures, where LVL is the number of levels in multigrid.  For each
%   level lvl, MG_DATA{lvl} contains at least the following fields :
%
%     error_ctrl : name of error functional to use to control the behavior
%         of the multigrid algorithm, see the ctrl option.
%
%     error_ctrl_flag : flag for error_ctrl nonempty.
%
%     error_ctrl_scale : see ctrl_scale option.
%
%     error_rel : flag for relative errors, see rel options.
%
%     error : list of error functionals.
%         error is a structure containing the error functional function
%         handles as fields.
%
%   Also, the options are stored in MG_DATA{1}.opt.error.
%
%   Examples :
%
%   1)  edit options directly
%
%       opt = mg_error;
%       opt.l2 = true;
%       opt.exact = true;
%       opt.ctrl = 'l2_iter';
%       opt.rel = true;
%       mg_data = mg_error(mg_data,opt);
%
%   2)  or set individual parameters
%
%       mg_data = mg_error(mg_data,...
%           'l2',true,'exact',true,'ctrl','l2_iter','rel',true);
%
%   3)  change individual parameters (after (1) or (2))
%
%       mg_data = mg_error(mg_data,'ctrl','energy_exact');
%
%   See also mg_mesh, mg_stima, mg_smooth, mg.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % return default options if no arguments are given
  if(nargin==0)
    error_opt.iter = true;
    error_opt.exact = false;
    
    error_opt.eucl = false;
    error_opt.l2 = false;
    error_opt.h1 = false;
    error_opt.energy = true;
    
    error_opt.custom = {};
    error_opt.custom_names = {};
    
    error_opt.ctrl = 'auto';
    error_opt.ctrl_scale = 'auto';
    error_opt.rel = false;
    
    error_opt.add_scapro = true;
    
    varargout{1} = error_opt;
    
    return
  end
  
  % read options
  if(length(varargin)==1)
    error_opt = varargin{1};
  else
    if(isfield(mg_data{1},'opt') && isfield(mg_data{1}.opt,'error'))
      error_opt = mg_data{1}.opt.error;
    else
      error_opt = mg_error;
    end
    for i=1:floor(0.5*length(varargin))
      error_opt.(varargin{2*i-1}) = varargin{2*i};
    end
  end
  
  % return options as second argument
  if(nargout>1)
    varargout{2} = error_opt;
  end
  
  % add options to first level of mg_data
  mg_data{1}.opt.error = error_opt;
  
  % determine assembly routines ...
  numvert = size(mg_data{1}.mesh.Elements,2);
  
  if(numvert==3) % triangular elements
    
    % ... for L2 norm
    if(error_opt.l2)
      assemL2 = @(mesh) assemMat_LFE(mesh,@MASS_LFE);
    end
    
    % ... for H1_0 norm
    if(error_opt.h1)
      assemH1 = @(mesh) assemMat_LFE(mesh,@STIMA_Lapl_LFE);
    end
    
  elseif(numvert==4) % quadrilateral elements
    
    % quadrature rule
    if(error_opt.l2 || error_opt.h1)
      quad = TProd(gauleg(0,1,2));
    end
    
    % ... for L2 norm
    if(error_opt.l2)
      assemL2 = @(mesh) assemMat_BFE(mesh,@MASS_BFE,quad);
    end
    
    % ... for H1_0 norm
    if(error_opt.h1)
      assemH1 = @(mesh) assemMat_BFE(mesh,@STIMA_Lapl_BFE,quad);
    end
  
  else
    error('Unknown mesh type');
  end
  
  % define names of error estimators
  names.eucl_iter = 'eucl_iter';
  names.l2_iter = 'l2_iter';
  names.h1_iter = 'h1_iter';
  names.energy_iter = 'energy_iter';
  
  names.eucl_exact = 'eucl_exact';
  names.l2_exact = 'l2_exact';
  names.h1_exact = 'h1_exact';
  names.energy_exact = 'energy_exact';
  
  % determine control error estimator
  ctrl_scale = 4;
  if(isequal(error_opt.ctrl,'auto'))
    name_ctrl = '';
    if(~isempty(error_opt.custom))
      if(~isempty(error_opt.custom_names))
        name_ctrl = error_opt.custom_names{1};
      else
        name_ctrl = 'custom1';
      end
    end
    if(error_opt.iter)
      if(error_opt.eucl)
        name_ctrl = names.eucl_iter;
      elseif(error_opt.l2)
        name_ctrl = names.l2_iter;
      elseif(error_opt.h1)
        name_ctrl = names.h1_iter;
        ctrl_scale = 2;
      elseif(error_opt.energy)
        name_ctrl = names.energy_iter;
        ctrl_scale = 2;
      end
    elseif(error_opt.exact)
      if(error_opt.eucl)
        name_ctrl = names.eucl_exact;
      elseif(error_opt.l2)
        name_ctrl = names.l2_exact;
      elseif(error_opt.h1)
        name_ctrl = names.h1_exact;
        ctrl_scale = 2;
      elseif(error_opt.energy)
        name_ctrl = names.energy_exact;
        ctrl_scale = 2;
      end
    end
  elseif(isa(error_opt.ctrl,'char'))
    name_ctrl = error_opt.ctrl;
  else
    error('Unknown control error estimator.');
  end 
  
  % override scaling behaviour of control error estimator
  if(isa(error_opt.ctrl_scale,'numeric'))
    ctrl_scale = error_opt.ctrl_scale;
  end
  
  % define error estimators on all levels
  for lvl=1:length(mg_data)
    
    % delete old error estimators
    if(isfield(mg_data{lvl},'error'))
      mg_data{lvl} = rmfield(mg_data{lvl},'error');
    end
    
    % initialize field error
    mg_data{lvl}.error = struct();
    
    % assemble matrices for inner products
    if(error_opt.iter || error_opt.exact)
      if(error_opt.l2)
        L2 = assemL2(mg_data{lvl}.mesh);
        L2 = L2(mg_data{lvl}.dofs,mg_data{lvl}.dofs);
        if(error_opt.add_scapro)
          mg_data{lvl}.L2 = L2;
        end
      end
      if(error_opt.h1)
        H1 = assemH1(mg_data{lvl}.mesh);
        H1 = H1(mg_data{lvl}.dofs,mg_data{lvl}.dofs);
        if(error_opt.add_scapro)
          mg_data{lvl}.H1 = H1;
        end
      end
    end
    
    % iteration errors
    if(error_opt.iter)
      if(error_opt.eucl)
        mg_data{lvl}.error.(names.eucl_iter) = @(x,y,varargin) sqrt(sum((x-y).^2));
      end
      if(error_opt.l2)
        mg_data{lvl}.error.(names.l2_iter) = @(x,y,varargin) sqrt((x-y)'*L2*(x-y));
      end
      if(error_opt.h1)
        mg_data{lvl}.error.(names.h1_iter) = @(x,y,varargin) sqrt((x-y)'*H1*(x-y));
      end
      if(error_opt.energy)
        mg_data{lvl}.error.(names.energy_iter) = @(x,y,varargin) sqrt((x-y)'*mg_data{lvl}.A*(x-y));
      end
    end
    
    % exact errors
    if(error_opt.exact)
      if(isempty(mg_data{lvl}.A))
        u = zeros(0,1);
      else
        u = mg_data{lvl}.A\mg_data{lvl}.b;
      end
      
      if(error_opt.eucl)
        mg_data{lvl}.error.(names.eucl_exact) = @(x,varargin) sqrt(sum((x-u).^2));
      end
      if(error_opt.l2)
        mg_data{lvl}.error.(names.l2_exact) = @(x,varargin) sqrt((x-u)'*L2*(x-u));
      end
      if(error_opt.h1)
        mg_data{lvl}.error.(names.h1_exact) = @(x,varargin) sqrt((x-u)'*H1*(x-u));
      end
      if(error_opt.energy)
        mg_data{lvl}.error.(names.energy_exact) = @(x,varargin) sqrt((x-u)'*mg_data{lvl}.A*(x-u));
      end
    end
    
    % custom error estimators
    for i=1:min(numel(error_opt.custom),numel(error_opt.custom_names))
      mg_data{lvl}.error.(error_opt.custom_names{i}) = error_opt.custom{i};
    end
    for j=i+1:numel(error_opt.custom)
      mg_data{lvl}.error.(sprintf('custom%.0f',j)) = error_opt.custom{j};
    end
    
    % error estimator to control multigrid
    mg_data{lvl}.error_ctrl = name_ctrl;
    mg_data{lvl}.error_ctrl_flag = ~isempty(name_ctrl);
    mg_data{lvl}.error_ctrl_scale = ctrl_scale;
    
    % check if this is a valid error functional
    if(~isfield(mg_data{lvl}.error,name_ctrl)...
        && ~isempty(fieldnames(mg_data{lvl}.error))...
        && ~isempty(name_ctrl))
      warning('LehrFEM:InvalidMGErrorFn',...
        '''%s'' is not a valid error functional on level %.0f.',name_ctrl,lvl);
    end
    
    % flag relative or full error
    mg_data{lvl}.error_rel = error_opt.rel;
    
  end % for loop over levels, index lvl
  
  % return mg_data multigrid data structure as first argument
  varargout{1} = mg_data;
  
return