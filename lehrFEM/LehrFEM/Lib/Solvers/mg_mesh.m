function varargout = mg_mesh(varargin)
%MG_MESH generate meshes for multigrid
%
%   MG_DATA = MG_MESH(OPT) generates the multigrid data structure MG_DATA
%   from the options OPT.  MG_DATA contains meshes with various amounts of
%   refinement and intergrid transfer operators.
%
%   OPT = MG_MESH() returns the default options wrapped in a structure.
%
%   [MG_DATA,OPT] = MG_MESH(OPT) returns both structures.
%
%   The options can be given in one of two forms: either as a data
%   structure similar to the return argument OPT, or by specifying
%   parameters and values, see the examples below.
%
%   The option parameters are
%
%     mesh : initial mesh data structure,
%         [default : unit square divided into two triangles].
%         The mesh may be composed of either triangles or quadrilaterals.
%
%     ref : 1-by-2 vector specifying the number of refinements of the
%         initial mesh to use for multigrid, [default : [3 7] ].
%       ref(1) : the coarsest mesh used by multigrid is constructed by
%         ref(1) regular refinements of the initial mesh.
%       ref(2) : the finest mesh used by multigrid is constructed by ref(2)
%         regular refinements of the initial mesh.
%
%     jiggle : non-negative double specifying the amount of jiggling in the
%         mesh refinements, [default : 0].  A positive value of about 1
%         leads to unstructured meshes that allow exact prolongation.
%
%     dir_flag : flag for Dirichlet boundary edges, [default : -1].
%
%     full : add full prolongation matrix to MG_DATA, [default : true].
%         The prolongation matrix MG_DATA{lvl}.P is restricted to the
%         non-Dirichlet-boundary vertices of the mesh.  If this argument is
%         true, then the full prologation matrix is stored in
%         MG_DATA{lvl}.P_full.
%
%     type : type of finite elements, [default : 'auto'].
%       'auto' : use linear or bilinear finite elements, depending on the
%         number of vertices in the elements of the initial mesh.
%       'LFE' : use linear finite elements.
%       'BFE' : use bilinear finite elements.
%       function handle : the function handle is used to contruct the
%         prolongation matrices; it must take the coarser mesh es the first
%         argument and the finer mesh as the second argument.
%
%   The multigrid data structure MG_DATA is a 1-by-LVL cell array of
%   structures, where LVL is the number of levels in multigrid,
%   ref(2)-ref(1)+1.  MG_DATA{1} corresponds to the coarsest level,
%   MG_DATA{LVL} to the finest.  For any level lvl, MG_DATA{lvl} contains
%   at least the following fields :
%
%     mesh : mesh data structure.
%
%     dofs : logical array specifying the non-Dirichlet-boundary degrees of
%         freedom.
%
%     n.all : the total number of degrees of freedom.
%
%     n.free : the number of non-Dirichlet-boundary degrees of freedom.
%
%   For lvl>1, the MG_DATA{lvl} alse contains
%
%     P : prolongation matrix from level lvl-1 to level lvl, for functions
%         that vanish on the Dirichlet boundary.
%
%   Also, the options are stored in MG_DATA{1}.opt.mesh.
%
%   Examples : 
%
%   1)  generate multigrid structure on L-shaped domain
%
%       Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%       opt = mg_mesh;
%       opt.mesh = Mesh;
%       opt.ref = [2 6];
%       mg_data = mg_mesh(opt);
%
%   2)  exactly the same thing can be achieved by
%
%       Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%       mg_data = mg_mesh('mesh',Mesh,'ref',[2 6]);
%
%   See also mg_stima, mg_smooth, mg_error, mg.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % return default options if no arguments are given
  if(nargin==0)
    mesh_opt.mesh.Coordinates = [0 0; ...
                                 1 0; ...
                                 1 1; ... 
                                 0 1];
    mesh_opt.mesh.Elements = [1 2 3; ...
                              1 3 4];
    mesh_opt.ref = [3 7];
    mesh_opt.jiggle = 0;
    mesh_opt.dir_flag = -1;
    mesh_opt.full = true;
    mesh_opt.type = 'auto';
    
    varargout{1} = mesh_opt;
    
    return
  end
    
  % read options
  if(length(varargin)==1)
    mesh_opt = varargin{1};
  else
    mesh_opt = mg_mesh;
    for i=1:floor(0.5*length(varargin))
      mesh_opt.(varargin{2*i-1}) = varargin{2*i};
    end
  end
  
  % return options as second argument
  if(nargout>1)
    varargout{2} = mesh_opt;
  end
  
  % initialize mesh
  CMesh = mesh_opt.mesh;
  
  CMesh = add_Edges(CMesh);
  CMesh = add_Edge2Elem(CMesh);
  CMesh.ElemFlag = zeros(size(CMesh.Elements,1),1);
  if(~isfield(CMesh,'BdFlags'))
    Loc = get_BdEdges(CMesh);
    CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
    CMesh.BdFlags(Loc) = mesh_opt.dir_flag;
  end
  
  for k=1:mesh_opt.ref(1)
    CMesh = refine_REG_jiggle(CMesh,mesh_opt.jiggle*0.1*2^-(mesh_opt.ref(1)-1));
  end
  
  % determine appropriate function for prolongation
  if(isa(mesh_opt.type,'function_handle'))
    getP = mesh_opt.type;
  else
    switch mesh_opt.type
      case 'LFE'
        getP = @get_PMat_LFE;
      case 'BFE'
        getP = @get_PMat_BFE;
      case 'auto'
        switch size(CMesh.Elements,2)
          case 3
            getP = @get_PMat_LFE;
          case 4
            getP = @get_PMat_BFE;
          otherwise
            error('Unknown mesh type');
        end
      otherwise
        error('Unknown prolongation type');
    end
  end
  
  % initialize mg_data multigrid data structure
  LVL = diff(mesh_opt.ref)+1;
  mg_data = cell(1,LVL);
  
  mg_data{1}.mesh = CMesh;
  
  mg_data{1}.dofs = true(size(CMesh.Coordinates,1),1);
  mg_data{1}.dofs(CMesh.Edges(CMesh.BdFlags==mesh_opt.dir_flag,:)) = false;
  
  mg_data{1}.n.all = size(CMesh.Coordinates,1);
  mg_data{1}.n.free = nnz(mg_data{1}.dofs);
  
  % construct finer meshes
  for lvl=2:LVL
    
    FMesh = refine_REG_jiggle(CMesh,mesh_opt.jiggle*0.1*2^-(LVL-lvl));
    mg_data{lvl}.mesh = add_MLevel(FMesh);
    mg_data{lvl}.n.all = size(FMesh.Coordinates,1);
    
    mg_data{lvl}.dofs = true(mg_data{lvl}.n.all,1);
    mg_data{lvl}.dofs(FMesh.Edges(FMesh.BdFlags==mesh_opt.dir_flag,:)) = false;
    mg_data{lvl}.n.free = nnz(mg_data{lvl}.dofs);
    
    P = getP(CMesh,FMesh);
    mg_data{lvl}.P = P(mg_data{lvl}.dofs,mg_data{lvl-1}.dofs);
    if(mesh_opt.full)
      mg_data{lvl}.P_full = P;
    end
    
    CMesh = FMesh;
  end
  
  % add options to all levels of mg_data
  mg_data{lvl}.opt.mesh = mesh_opt;
  
  % return mg_data multigrid data structure as first argument
  varargout{1} = mg_data;
  
return