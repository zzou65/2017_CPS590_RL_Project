% Generates a unstructured mesh of a nose-shaped domain

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  X0 = [0.875 -0.49607837];                                                         % Intersection point of circles
  X1 = [0.875  0.49607837];                                                         % Intersection point of circles                                          
  BBOX = [-1 -1; 2 1];                                                              % Bounding box
  H0 = 0.1;                                                                        % Initial mesh width
  DHANDLE = inline('dist_union(dist_circ(x,[0 0],1),dist_circ(x,[1 0],0.5))','x');  % Signed distance function
  HHANDLE = @h_uniform;                                                             % Element size function
  FIXEDPOS = [X0; X1];                                                              % Fixed boundary vertices of the mesh
  DISP = 1;                                                                         % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);  
  Mesh = add_Edges(Mesh);
  
  % Generate plots
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;