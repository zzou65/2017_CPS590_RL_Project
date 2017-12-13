% Generates figures of meshes of the unit circle:
%  MeshDiskd1.eps,  MeshDiskd2.eps,

% Copyright 2005-2005 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants

  DHANDLE = @dist_circ;  % Signed distance function for a circle
  C = [0 0];             % Center of circle
  R = 1;                 % Radius of circle
  INIT_NREFS = 1;        % Number of initial mesh refinements
  NREFS = 2;             % Total number mesh refinements

  % Initialize mesh

  Mesh = load_Mesh('Coord_Ball.dat','Elem_Ball.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;

  % Do initial mesh refinements

  for i = 1:INIT_NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
  end

  % Generate .eps files

  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHANDLE,C,R);
    fig = plot_Mesh(Mesh,'as');
%     filename = ['MeshDiskd' int2str(i) '.eps'];
% %     print('-depsc',filename);
%     close(fig);
%     system(['gv ' filename ' &']);
  end

  % Clear memory

  clear all;
