% Generates a triangular unstrucutred mesh of a star with a circular hole

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  theta = linspace(0,2*pi,7);
  theta = theta(1:6)';

  % Initialize constants
  
  BBOX = [-1 -1; 1 1];                                                                                % Bounding box
  H0 = 0.05;                                                                                          % Initial mesh width
  DHANDLE = inline(['dist_diff(dist_union(dist_tri(x,[0 1],[-sqrt(3)/2 -1/2],[sqrt(3)/2 -1/2]),' ...  % Signed distance function
                                         'dist_tri(x,[0 -1],[sqrt(3)/2 1/2],[-sqrt(3)/2 1/2])),' ...
                                         'dist_circ(x,[0 0],.5))'],'x');
  HHANDLE = @h_uniform;                                                                                                                                                           % Element size function
  FIXEDPOS = [cos(theta)/sqrt(3) sin(theta)/sqrt(3); ...                                              % Fixed positions of the mesh                                                                                                                          
              cos(theta+pi/6)    sin(theta+pi/6)];                                        
  DISP = 1;                                                                                           % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
    
  % Generate plots of the mesh
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  % Clear memory
  
  clear all;
  