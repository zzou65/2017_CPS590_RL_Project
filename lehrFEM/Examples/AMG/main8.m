function [] = main8()
% coarse grids for regular refinements
%
%   Plots AMG coarse grid points for regularly refined quadrilateral and
%   triangular meshes.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  %%% Triangular mesh

  %% Initialize

  % Generate coarse mesh

  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % Refine mesh
  
  for i = 1:2
    Mesh = refine_REG(Mesh);
  end
  
  % Construct stiffness matrix (ignoring boundary)
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
  
  %% Construct AMG Data

  % Define Options
  
  opt = AMGDefaultOptions;
  opt.mincoarsegrid = 10;
  opt.levsmax = 0;
  
  % Generate multilevel structure
  
  data = AMGSetup(A,opt);
  
  %% plot coarse grid points
  
  plot_AMG_coarse(Mesh,data);
  title('\bf Coarse Grid Points');
  

  %%% Quadrilateral mesh

  %% Initialize

  % Generate coarse mesh

  Mesh = load_Mesh('Coord_Sqr_QElem.dat','Elem_Sqr_QElem.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % Refine mesh
  
  for i = 1:2
    Mesh = refine_REG(Mesh);
  end
  
  % Construct stiffness matrix (ignoring boundary)
  
  A = assemMat_BFE(Mesh,@STIMA_Lapl_BFE,TProd(gauleg(0,1,2)));
  
  
  %% Construct AMG Data

  % Define Options
  
  opt = AMGDefaultOptions;
  opt.mincoarsegrid = 10;
  opt.levsmax = 0;
  
  % Generate multilevel structure
  
  data = AMGSetup(A,opt);
  
  %% plot coarse grid points
  
  plot_AMG_coarse(Mesh,data);
%   title('\bf Coarse Grid Points');
