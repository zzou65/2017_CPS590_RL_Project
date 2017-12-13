% Generates figures of unstructured meshes of the unit square:
%  MeshSquare1.eps,  MeshSquare2.eps,

% Copyright 2005-2005 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants

  INIT_NREFS = 1;        % Number of initial mesh refinements
  NREFS = 2;             % Total number mesh refinements

  % Initialize mesh

  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;

  % Do initial mesh refinements

  for i = 1:INIT_NREFS
    Mesh = refine_REG(Mesh);
  end

  % Generate .eps files

  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
    Loc = get_BdEdges(Mesh);
    Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
    FixedPos = zeros(size(Mesh.Coordinates,1),1);
    FixedPos(Loc) = 1;
    Mesh = jiggle(Mesh,FixedPos);
    fig = plot_Mesh(Mesh,'as');
    filename = ['MeshSquare' int2str(i) '.eps'];
    print('-depsc',filename);
    close(fig);
    system(['gv ' filename ' &']);
  end
  
  % Clear memory
  
  clear all;
