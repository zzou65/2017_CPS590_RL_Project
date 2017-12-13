% Refine mesh at each corner node and plot out distribution of polynomial
% degrees.

% Copyright 2005-2006 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constant
    
  NREFS = 5;

  % Initialize mesh
    
  Mesh = load_Mesh('Coord.dat','Elem.dat');
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags(Loc) = -1;
    
  CNodes = (1:9)';
  
  % Adaptive mesh reefinement
  
  Mesh = init_LEB(Mesh);
  for i = 1:NREFS
    Mesh = refine_hp(Mesh,CNodes);  
  end
  
  Mesh_hp.Coordinates = Mesh.Coordinates;
  Mesh_hp.Elements = Mesh.Elements;
  Mesh_hp.ElemFlag = zeros(size(Mesh_hp.Elements,1),1);
  Mesh_hp = add_Edges(Mesh_hp);
  Loc = get_BdEdges(Mesh_hp);
  Mesh_hp.BdFlags = zeros(size(Mesh_hp.Edges,1),1);
  Mesh_hp.BdFlags(Loc) = -1;
  Mesh_hp = add_Edge2Elem(Mesh_hp);
  
  % Assign polynomial degrees to elements
  
  [EDofs,CDofs,ElemDeg] = assign_pdeg(Mesh_hp,CNodes,NREFS);
  
  % Assign number to elements

  fig = figure('Name','Polynomial degrees');  
  patch('Faces',Mesh.Elements, ...
        'Vertices',Mesh.Coordinates, ...
        'FaceVertexCData',ElemDeg, ...
        'EdgeColor','k', ...
        'FaceColor','flat');
  set(gca,'CLim',[1 NREFS],'DataAspectRatio',[1 1 1]);
  colormap(jet);
  alpha(.9);
  colorbar;
  set(gcf,'renderer','openGL');
  print('-depsc', 'refine_hp.eps');
    
  % Clear memory
    
  clear all;
    