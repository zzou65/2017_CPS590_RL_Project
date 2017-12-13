% Generates an triangular unstructured mesh of the cross section of a
% quadrilateral pipe

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  X0 = [-1 -1];
  XI = [-0.5 -0.5]% Lower left point the pipe
  A = 2;
  AI = 1;% Side length of the pipe
  BBOX = [-1 -1; 1 1];                                                                % Bounding box
  H0 = 0.25;                                                                          % Initial mesh width
  DHANDLE = inline('dist_diff(dist_rect(x,[-1 -1],2,2),dist_rect(x,[-0.5 -0.5],1,1))','x');  % Signed distance function
  HHANDLE = @h_uniform;                                                               % Element size function
  FIXEDPOS = [X0; X0+[A 0]; X0+[A A]; X0+[0 A];XI; XI+[AI 0]; XI+[AI AI]; XI+[0 AI]];     % Fixed boundary vertices of the mesh  
  DISP = 1;                                                                           % Display flag
  
  % Generate mesh
  
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  Mesh = add_Edges(Mesh);
  
  % Generate plots
  
  plot_Mesh(Mesh,'as');
  plot_Qual(Mesh);
  plot_USR(Mesh);
  
  save_Mesh(Mesh,'PMLSquareCoord.dat','PMLSquareElem.dat');
  
  % Clear memory
  
  clear all;