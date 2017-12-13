% Generates an triangular unstructured mesh of the unit square

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  X0 = [-1 -1];                                   % Lower left corner point of rectangle
  A = 2;                                          % Length of rectangle
  BBOX = [X0; X0+[A A]];                          % Bounding box
  H0 = 0.2;                                       % Initial mesh width
  DHANDLE = @dist_rect;                           % Signed distance function
  HHANDLE = @h_uniform;                           % Element size function
  FIXEDPOS = [X0; X0+[A 0]; X0+[A A]; X0+[0 A]];  % Fixed boundary vertices of the mesh  
  DISP = 1;                                       % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,X0,A,A);
  Mesh = add_Edges(Mesh);
  
  % Generate plots
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
 % clear all;