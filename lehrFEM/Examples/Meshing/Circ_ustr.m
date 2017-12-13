% Generates a triangular unstrucutured mesh of the unit circle

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  C = [0 0];             % Center of circle
  R = 1;                 % Radius of circle 
  BBOX = [-1 -1; 1 1];   % Bounding box
  H0 = 0.1;              % Initial mesh width
  DHANDLE = @dist_circ;  % Signed distance function
  HHANDLE = @h_uniform;  % Element size function
  FIXEDPOS = [];         % Fixed boundary vertices of the mesh
  DISP = 1;              % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);  
  Mesh = add_Edges(Mesh);
  
  % Generate plots
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;
  