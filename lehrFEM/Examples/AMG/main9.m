function [] = main9()
% coarse grid points for local refinements

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % define AMG options
  
  opt = AMGDefaultOptions;
  opt.mincoarsegrid = 10;
  opt.levsmax = 0;

  %% Refinement towards a point

  % initialize mesh

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:0
    Mesh = refine_REG(Mesh);
  end
  
  % generate refined mesh using local multigrid
  
  [U,rMesh,flag,iter,ndof,err,timer,dofslvl,A] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,@pointref,0.5,1,1,1,@gs_smooth,0,8);
  
  % generate multilevel structure
  
  data = AMGSetup(A,opt);
  
  % plot coarse grid points
  
  plot_AMG_coarse(rMesh,data);
  
  
  %% Refinement towards a circle
  
  % initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr_Big.dat','Elem_Sqr_Big.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:0
    Mesh = refine_REG(Mesh);
  end
  
 
  % generate refined mesh using local multigrid
  
  assemLoad = @(mesh) assemLoad_LFE_S1(mesh,gauleg(0,1,2));
  g_D = @(x,varargin) zeros(size(x,1),1);

  [U,rMesh,flag,iter,ndof,err,timer,dofslvl,A] = locmg_solve([],Mesh,{assemLoad},g_D,@circref,0.5,1,1,[1,1],@gs_smooth,0.05,8);

  % generate multilevel structure
  
  data = AMGSetup(A,opt);
  
  % plot coarse grid points
  
  plot_AMG_coarse(rMesh,data);

return