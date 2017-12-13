% Morphs a triangular mesh of the difference of two circles into a
% quadrilateral mesh
  
%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  BBOX = [0 0; 1 1];                                                                                                  % Bounding box
  H0 = 0.1;                                                                                                           % Initial mesh width
  DHANDLE = inline('dist_diff(dist_diff(dist_diff(dist_circ(x,[0 0],1),x(:,2)),x(:,1)),dist_circ(x,[0 0],.5))','x');  % Signed distance function
  HHANDLE = @h_uniform;                                                                                               % Element size function
  FIXEDPOS = [0.5 0; 0 0.5; 1 0; 0 1];                                                                                % Fixed boundary vertices of the mesh
  DISP = 1;                                                                                                           % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP); 
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
  