% Morphs a triangular mesh of the L-shaped domain into a quadrilateral mesh

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 3;                                                                            % Number of unifrom red refinements
  DHANDLE = inline('dist_diff(dist_rect(x,[-1 -1],2,2),dist_rect(x,[0 -1],1,1))','x');  % Signed distance function
  
  % Load mesh from file

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  
  % Add edge data structure
  
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do NREFS uniform refinement steps
  
  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHANDLE);
  end
  plot_Qual(Mesh);
  
  % Morph the mesh
  
  Mesh = morph(Mesh);
  
  % Laplacian smoothing
  
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
  FixedPos = zeros(size(Mesh.Coordinates,1),1);
  FixedPos(Loc) = 1;
  Mesh = smooth(Mesh,FixedPos);
  plot_Qual(Mesh);
  
  % Clear memory
  
  clear all;
  
  