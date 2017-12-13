% Generates a triangular unstrucutured mesh of a dumbbell

%   Copyright 2005-2005 Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  BBOX = [-3 -1;3 1];                                                                                                                % Bounding box
  H0 = 0.2;                                                                                                                          % Initial mesh width
  DHANDLE = inline('dist_union(dist_union(dist_rect(x,[-3 -1],2.5,2),dist_rect(x,[.5 -1],2.5,2)),dist_rect(x,[-1 -.2],2,.4))','x');  % Signed distance 
  HHANDLE = @h_uniform;                                                                                                              % Element size function
  FIXEDPOS = [-3 -1;-.5 -1;-3 1; -.5 1;-.5 .2;-.5 -.2;.5 -.2;.5 .2;.5 -1;.5 1;3 1;3 -1;];                                            % Fixed boundary vertices
  DISP = 1;                                                                                                                          % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  Mesh = add_Edges(Mesh);
  
  % Generate figure
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;