% Generates a quadrilateral structured mesh of the unit square

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
 
  % Initialize constants
  
  NREFS = 3;  % Number of unifrom red refinements
 
  % Load mesh from file
   
  Mesh = load_Mesh('Coord_Sqr_QElem.dat','Elem_Sqr_QElem.dat');
  
  % Add edge data structure
  
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  
  % Do NREFS uniform refinement steps
   
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end    
  plot_Mesh(Mesh,'setpa');
  
  % Clear memory
  
  clear all;
  