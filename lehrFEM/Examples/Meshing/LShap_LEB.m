% Generates a mesh of the L-shaped domain by largest edge bisection from
% an initial coarse mesh

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  NREFS = 15;

  % Load mesh from file

  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  
  % Add edge data structure
  
  Mesh.ElemFlag = [1 2 2 1]';
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
 
  % Initialize largest edge bisection
  
  Mesh = init_LEB(Mesh);  
  
  % Do largest edge bisection
  
  for i = 1:NREFS
    Mesh = refine_LEB(Mesh,[1 2 3 4 5]);
  end  
  plot_Mesh(Mesh,'as');
  
  % Clear memory
  
  clear all;