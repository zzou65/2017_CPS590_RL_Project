% Generates figures of structured meshes of the L-shaped domain:
%  Lshape0.eps,      Lshape1.eps,      Lshape2.eps,
%  MeshLshape1.eps,  MeshLshape2.eps.

% Copyright 2005-2005 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  INIT_NREFS = 1;  % Number of initial mesh refinements
  NREFS = 3;       % Total number mesh refinements
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
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
    fig = plot_Mesh(Mesh,'as');
    filename = ['Lshape' int2str(i-1) '.eps'];
    print('-depsc',filename);
    system(['gv ' filename ' &']);
    if(i == 1 || i == 2)
      filename = ['MeshLshape' int2str(i) '.eps'];
      print('-depsc',filename);
      system(['gv ' filename ' &']);
    end
    close(fig);
  end
  
  % Clear memory
  
  clear all;
  