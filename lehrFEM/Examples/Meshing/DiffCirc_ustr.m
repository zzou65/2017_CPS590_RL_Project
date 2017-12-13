% Generates a triangular unstrucutured mesh of the difference of two circles

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  BBOX = [0 0; 1 1];                                                                                                  % Bounding box
  H0 = 0.05;                                                                                                           % Initial mesh width
  DHANDLE = inline('dist_diff(dist_diff(dist_diff(dist_circ(x,[0 0],1),x(:,2)),x(:,1)),dist_circ(x,[0 0],.5))','x');  % Signed distance function
  
  
  HHANDLE = @h_uniform;                                                                                               % Element size function
  FIXEDPOS = [0 .5;.5 0; 1 0; 0 1];                                                                                  % Fixed boundary vertices of the mesh
  DISP = 1;                                                                                                           % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  Mesh = add_Edges(Mesh);
  
  % Generate figure
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  %clear all;