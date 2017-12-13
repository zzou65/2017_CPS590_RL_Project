function varargout = mg_stima(mg_data,varargin)
%MG_STIMA assemble stiffness matrices and load vectors for multigrid
%   
%   MG_Data = MG_STIMA(MG_DATA,OPT) adds stiffness matrices and (possibly)
%   load vectors to the multigrid data structure MG_DATA.  MG_DATA should
%   contain mesh information.
%
%   OPT = MG_STIMA() returns the default options wrapped in a structure.
%
%   [MG_DATA,OPT] = MG_STIMA(MG_DATA,OPT) returns both structures.
%
%   The options can be given in one of two forms: either as a data
%   structure similar to the return argument OPT, or by specifying
%   parameters and values, see the examples below.  In the latter
%   case, if MG_DATA already contains stiffness matrices, it is possible to
%   include only the options one wishes to change.
%
%   The option parameters are
%
%     stima : type of assembly routine for the stiffness matrix, values :
%       'auto' [default] : use standard assembly routines for linear or
%         bilinear finite elements.
%       'elem' : used standard assembly routine with element contributions
%         given by the function handle stima_assem with parameters given by
%         the cell array stima_param.
%       'assem' : use the assembly routine given by the funtion handle
%         stima_assem; stima_assem must take the single argument mesh.
%
%     stima_assem : alternative routine for assebling the stiffness matrix,
%       see stima.
%
%     stima_param : parameters for alternative assembly routine, see stima.
%
%     stima_quad : quadrature rule to use for stiffness matrix, values:
%       'auto' [default] : use P7O6() for linear elements and 
%           TProd(gauleg(0,1,2)) for bilinear elements.
%       quadrule : use quadrature rule quadrule.
%
%     add_load : add load vectors to MG_DATA, [default : true].
%
%     load : assembly routine for load vector, values : 
%       'auto' [default] : use standard assembly routines.
%       (function handle) : use load as the assembly routine.  It must take
%         the single argument mesh.
%
%     load_quad : quadrature rule to use for load vector, values :
%       'auto' [default] : use P7O6 for LFE and TProd(gauleg(0,1,2)) for
%         BFE.
%       quadrule : use quadrature rule quadrule.
%
%     full : add full stiffness matrices and load vectors to data
%       structure, [default : false].
%
%     dir_flag : flag for Dirichlet boundary edges, [default : -1].
%
%     f : function handle for right hand side,
%       [default : @(x,varargin) zeros(size(x,1),1)].
%
%     gd : function handle for Dirichlet boundary data,
%       [default : @(x,varargin) zeros(size(x,1),1)].
%
%     c : function handle or scalar for (scalar) conductivity,
%       [default : 1].
%       
%   The multigrid data structure MG_DATA is a 1-by-LVL cell array of
%   structures, where LVL is the number of levels in multigrid.  When
%   constructed with default arguments, MG_DATA{lvl} contains at least the 
%   following fields for any level lvl :
%
%     A : stiffness matrix (for free degrees of freedom).
%
%     u_bd : boundary data for the solution vector.
%
%     b : load vector (for free degrees of freedom).
%
%   MG_DATA{1}, which corresponds to the coarsest level, also  contains
%
%     L, U : LU-decomposition of A, ie. A = L*U.
%
%   Also, the options are stored in MG_DATA{1}.opt.stima.
%
%   Examples : 
%
%   1)  edit options directly
%
%       opt = mg_stima;
%       opt.f = @(x,varargin) -4*ones(size(x,1),1);
%       opt.gd = @(x,varargin) x(:,1).^2+x(:,2).^2;
%       mg_data = mg_stima(mg_data,opt);
%
%   2)  or set individual parameters
%
%       mg_data = mg_stima(mg_data,...
%         'f',@(x,varargin) -4*ones(size(x,1),1),...
%         'gd',@(x,varargin) x(:,1).^2+x(:,2).^2);
%
%   3)  change individual parameters (after (1) or (2))
%
%       mg_data = mg_stima(mg_data,'gd',@(x,varargin) ones(size(x,1),1));
%
%   See also mg_mesh, mg_smooth, mg_error, mg.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % return default options if no arguments are given
  if(nargin==0)
    stima_opt.stima = 'auto';
    stima_opt.stima_param = {};
    stima_opt.stima_assem = [];
    stima_opt.stima_quad = 'auto';
    
    stima_opt.add_load = true;
    stima_opt.load = 'auto';
    stima_opt.load_quad = 'auto';
    
    stima_opt.full = false;
    stima_opt.dir_flag = -1;
    
    stima_opt.f = @(x,varargin) zeros(size(x,1),1);
    stima_opt.gd = @(x,varargin) zeros(size(x,1),1);
    stima_opt.c = 1;
    
    varargout{1} = stima_opt;
    
    return
  end
  
  % read options
  if(length(varargin)==1)
    stima_opt = varargin{1};
  else
    if(isfield(mg_data{1},'opt') && isfield(mg_data{1}.opt,'stima'))
      stima_opt = mg_data{1}.opt.stima;
    else
      stima_opt = mg_stima;
    end
    for i=1:floor(0.5*length(varargin))
      stima_opt.(varargin{2*i-1}) = varargin{2*i};
    end
  end
  
  % return options as second argument
  if(nargout>1)
    varargout{2} = stima_opt;
  end
  
  % add options to first level of mg_data
  mg_data{1}.opt.stima = stima_opt;
  
  % determine assembly routines
  numvert = size(mg_data{1}.mesh.Elements,2);
  
  if(numvert==3) % triangular elements
    
    % quadrature rules
    if(isequal(stima_opt.stima_quad,'auto'))
      stima_quad = P7O6();
    else
      stima_quad = stima_opt.stima_quad;
    end
    if(stima_opt.add_load)
      if(isequal(stima_opt.load_quad,'auto'))
        load_quad = P7O6();
      else
        load_quad = stima_opt.load_quad;
      end
    end
    
    % ... for stiffness matrix
    if(isequal(stima_opt.stima,'auto'))
      if(isa(stima_opt.c,'function_handle'))
        assemMat = @(mesh) assemMat_LFE(mesh,@STIMA_Heat_LFE,stima_quad,stima_opt.c);
      elseif(isa(stima_opt.c,'numeric'))
        assemMat = @(mesh) stima_opt.c*assemMat_LFE(mesh,@STIMA_Lapl_LFE);
      end
    elseif(isequal(stima_opt.stima,'elem'))
      assemMat = @(mesh) assemMat_LFE(mesh,stime_opt.stima_assem,stime_opt.stima_param{:});
    elseif(isequal(stima_opt.stima,'assem'))
      assemMat = stima_opt.stima_assem;
    else
      error('Unknown assembly routine for stiffness matrix.');
    end
    
    % ... for load vector
    if(stima_opt.add_load)
      if(isequal(stima_opt.load,'auto'))
        assemLoad = @(mesh) assemLoad_LFE(mesh,load_quad,stima_opt.f);
      elseif(isa(stima_opt.load,'function_handle'))
        assemLoad = stima_opt.load;
      else
        error('Unknown assembly routine for load vector.');
      end
    end
    
    % ... for Dirichlet boundary conditions
    assemDir = @(mesh) assemDir_LFE(mesh,stima_opt.dir_flag,stima_opt.gd);
      
  elseif(numvert==4) % quadrilateral elements
    
    % quadrature rules
    if(isequal(stima_opt.stima_quad,'auto'))
      stima_quad = TProd(gauleg(0,1,2));
    else
      stima_quad = stima_opt.stima_quad;
    end
    if(stima_opt.add_load)
      if(isequal(stima_opt.load_quad,'auto'))
        load_quad = TProd(gauleg(0,1,2));
      else
        load_quad = stima_opt.load_quad;
      end
    end
    
    % ... for stiffness matrix
    if(isequal(stima_opt.stima,'auto'))
      if(isa(stima_opt.c,'function_handle'))
        assemMat = @(mesh) assemMat_BFE(mesh,@STIMA_Heat_BFE,stima_quad,stima_opt.c);
      elseif(isa(stima_opt.c,'numeric'))
        assemMat = @(mesh) stima_opt.c*assemMat_BFE(mesh,@STIMA_Lapl_BFE,stima_quad);
      end
    elseif(isequal(stima_opt.stima,'elem'))
      assemMat = @(mesh) assemMat_BFE(mesh,stime_opt.stima_assem,stime_opt.stima_param{:});
    elseif(isequal(stima_opt.stima,'assem'))
      assemMat = stima_opt.stima_assem;
    else
      error('Unknown assembly routine for stiffness matrix.');
    end
    
    % ... for load vector
    if(stima_opt.add_load)
      if(isequal(stima_opt.load,'auto'))
        assemLoad = @(mesh) assemLoad_BFE(mesh,load_quad,stima_opt.f);
      elseif(isa(stima_opt.load,'function_handle'))
        assemLoad = stima_opt.load;
      else
        error('Unknown assembly routine for load vector.');
      end
    end
    
    % ... for Dirichlet boundary conditions
    assemDir = @(mesh) assemDir_BFE(mesh,stima_opt.dir_flag,stima_opt.gd);
    
  else
    error('Unknown mesh type');
  end
  
  % assemble stiffness matrix and load vector on coarsest grid
  A = assemMat(mg_data{1}.mesh);
  mg_data{1}.A = A(mg_data{1}.dofs,mg_data{1}.dofs);
  if(stima_opt.full)
    mg_data{1}.A_full = A;
  end
  mg_data{1}.u_bd = assemDir(mg_data{1}.mesh);
  if(stima_opt.add_load)
    b = assemLoad(mg_data{1}.mesh);
    b_dir = b - A*mg_data{1}.u_bd;
    mg_data{1}.b = b_dir(mg_data{1}.dofs);
    if(stima_opt.full)
      mg_data{1}.b_full = b;
    end
  end
  
  % calculate LU-factorization of A on the coarsest level
  [mg_data{1}.L,mg_data{1}.U] = lu(mg_data{1}.A);
  
  % assemble stiffness matrices and load vectors on finer levels
  for lvl=2:length(mg_data)
    A = assemMat(mg_data{lvl}.mesh);
    mg_data{lvl}.A = A(mg_data{lvl}.dofs,mg_data{lvl}.dofs);
    if(stima_opt.full)
      mg_data{lvl}.A_full = A;
    end
    mg_data{lvl}.u_bd = assemDir(mg_data{lvl}.mesh);
    if(stima_opt.add_load)
      b = assemLoad(mg_data{lvl}.mesh);
      b_dir = b - A*mg_data{lvl}.u_bd;
      mg_data{lvl}.b = b_dir(mg_data{lvl}.dofs);
      if(stima_opt.full)
        mg_data{lvl}.b_full = b;
      end
    end
  end
  
  % return mg_data multigrid data structure as first argument
  varargout{1} = mg_data;
  
return