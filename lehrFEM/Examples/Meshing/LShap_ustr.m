% Generates a triangular unstrucutred mesh of the L-shaped domain

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  BBOX = [-1 -1; 1 1];                                                                  % Bounding box
  H0 = 0.1;                                                                             % Initial mesh width
  DHANDLE = inline('dist_diff(dist_rect(x,[-1 -1],2,2),dist_rect(x,[0 -1],1,1))','x');  % Signed distance function
  HHANDLE = @h_uniform;                                                                 % Element size function
  FIXEDPOS = [-1 -1; 0 -1; 0 0; 1 0; 1 1; -1 1];                                        % Fixed boundary vertices of the mesh
  DISP = 1;                                                                             % Display flag
  
  % Generate mesh
  
  disp('Generating mesh for L-shaped domain');
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  % plot_Mesh(Mesh,'as');
  disp('Adding edges');
  Mesh = add_Edges(Mesh);
  
  % Generate plots of the mesh
  
  disp('Plotting mesh with edges');
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;
  