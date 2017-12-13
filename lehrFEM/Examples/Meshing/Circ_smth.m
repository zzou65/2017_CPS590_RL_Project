% Smoothes a structured mesh of the unit circle

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 3;             % Number of unifrom red refinements
  DHANDLE = @dist_circ;  % Signed distance function
  C = [0 0];             % Center of the circle
  R = 1;                 % Radius of the circle
  
  % Load mesh from file

  Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
  Mesh = shift(Mesh,C);
  
  % Add edge data structure
  
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do NREFS uniform refinement steps
   
  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
  end  
  plot_Qual(Mesh);
  
  % Compute fixed positions of the mesh

  Loc = get_BdEdges(Mesh);
  FixedPos = zeros(size(Mesh.Edges,1),1);
  FixedPos(unique([Mesh.Edges(Loc,1) Mesh.Edges(Loc,2)])) = 1;
  
  % Laplacian smoothing
  
  Mesh = smooth(Mesh,FixedPos);
  plot_Qual(Mesh);
  
  % Clear memory
  
  clear all;
  