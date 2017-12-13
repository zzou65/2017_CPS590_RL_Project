% Generates figures of unstructured meshes of a pipe:
%  MeshPipe1.eps,  MeshPipe2.eps,

% Copyright 2005-2005 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants

  H0 = 0.25;   % Initial mesh width
  NREFS = 2;   % Total number mesh refinements
  R_OUT = 1;   % Radius of outer circle
  R_IN = 0.5;  % Radius of inner circle                                                                                    
  
  % Initialize mesh
  
  DHandle = inline(['dist_diff(dist_circ(x,[0 0],' num2str(R_OUT) ...
                    '),dist_circ(x,[0 0],' num2str(R_IN) '))'],'x');
  Mesh = init_Mesh([-R_OUT -R_OUT; R_OUT R_OUT],H0,DHandle,@h_uniform,[],1);
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);  
  Mesh.BdFlags(Loc) = -1;

  % Generate .eps files

  for i = 1:NREFS
    Mesh = refine_REG(Mesh,DHandle);
    fig = plot_Mesh(Mesh,'as');
    filename = ['MeshPipe' int2str(i) '.eps'];
    print('-depsc',filename);
    close(fig);
    system(['gv ' filename ' &']);
  end
  
  % Clear memory
  
  clear all;
