% Generates a triangular unstrucutred mesh of the L-shaped domain

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  BBOX = [-2 -2; 2 2];                                       % Bounding box
  H0 = 0.1;                                                  % Initial mesh width
  DHANDLE = inline('dist_tri(x,[-2 2],[-1 -1],[2 0])','x');  % Signed distance function
  HHANDLE = @h_uniform;                                      % Element size function
  FIXEDPOS = [-2 2;-1 -1;2 0];                               % Fixed boundary vertices of the mesh
  DISP = 1;                                                  % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  plot_Mesh(Mesh,'as');
  Mesh = add_Edges(Mesh);
  
  % Generate plots of the mesh
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;
  