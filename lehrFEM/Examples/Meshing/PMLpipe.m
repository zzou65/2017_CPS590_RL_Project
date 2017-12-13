% Generates an triangular unstructured mesh of the cross section of a
% quadrilateral pipe

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  X0 = [-1 -1];% Lower left point the pipe
  X1 = [-2 -2];
  B  = 4;
  A = 2;                                                                              % Side length of the pipe
  BBOX = [-2 -2; 2 2];                                                                % Bounding box
  H0 = 0.15;                                                                          % Initial mesh width
  DHANDLE = inline('dist_diff(dist_rect(x,[-1 -1],2,2),dist_circ(x,[0 0],.5))','x');  % Signed distance function
  DHANDLE_PML = inline('dist_diff(dist_rect(x,[-2 -2],4,4), dist_rect(x,[-1 -1],2,2))','x');  % Signed distance function
  DHANDLE_PML1 = inline('dist_rect(x,[-1 -2],2,1)','x');  % Signed distance function
  DHANDLE_PML2 = inline('dist_rect(x,[-1 -1],1,2)','x');  % Signed distance function
  DHANDLE_PML3 = inline('dist_rect(x,[-1 1],2,1)','x');  % Signed distance function
  DHANDLE_PML4 = inline('dist_rect(x,[-2 -1],1,2)','x');  % Signed distance function
  DHANDLE_PML5 = inline('dist_rect(x,[-2 -2],1,1)','x');  % Signed distance function
  DHANDLE_PML6 = inline('dist_rect(x,[-2 1],1,1)','x');  % Signed distance function
  DHANDLE_PML7 = inline('dist_rect(x,[1 1],1,1)','x');  % Signed distance function
  DHANDLE_PML8 = inline('dist_rect(x,[-2 1],1,1)','x');  % Signed distance function
  DHANDLE_all = inline('dist_diff(dist_rect(x,[-2 -2],4,4),dist_circ(x,[0 0],.5))','x');  % Signed distance function
  HHANDLE = @h_uniform;                                                               % Element size function
  FIXEDPOS = [X0; X0+[A 0]; X0+[A A]; X0+[0 A]] ;% Fixed boundary vertices of the mesh  
  FIXEDPOS_PML = [X1; X1+[B 0]; X1+[B B]; X1+[0 B]] ;
  DISP = 1;                                                                           % Display flag
  NREFS=3;
    
  % Generate mesh
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Fixed=[];
  for i=1:size(Loc)
      if (norm(Mesh.Coordinates(Mesh.Edges(i,1)))>0.5)
         Fixed=[Fixed; i];
      end
  end
  % Generate PML-mesh
  FIXEDPOS_PML=[Mesh.Coordinates(unique(Mesh.Edges(Fixed,:)),:);FIXEDPOS_PML];
  PMLMesh = init_Mesh(BBOX,H0,DHANDLE_PML,HHANDLE,FIXEDPOS_PML, DISP);
  PMLMesh = add_Edges(PMLMesh);
  Loc = get_BdEdges(PMLMesh);
  PMLMesh.BdFlags = zeros(size(PMLMesh.Edges,1),1);
  PMLMesh.BdFlags(Loc) = -1;

  % Generate plots
  plot_Mesh(PMLMesh,'as');
  
  allMesh=merge_Mesh(Mesh,PMLMesh);
  Loc = get_BdEdges(allMesh);
  allMesh.BdFlags = zeros(size(allMesh.Edges,1),1);
  allMesh.BdFlags(Loc) = -1;
  plot_Mesh(allMesh);
  for i=1:NREFS
      Mesh=refine_REG(allMesh,DHANDLE_all);
      plot_Mesh(Mesh,'as');
  end
  
  % Clear memory
  
  clear all;
  