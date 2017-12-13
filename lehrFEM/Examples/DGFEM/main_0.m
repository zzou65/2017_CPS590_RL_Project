% Plot mesh inculding edge normals.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  clear all
  NREFS = 2;  % Number of red refinement steps
  
  % Initialize mesh
  
  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 3 4];
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  Mesh = rmfield(Mesh,'BdFlags');
  
  % Plot normals
  
  MidPts = (Mesh.Coordinates(Mesh.Edges(:,1),:) + ...
            Mesh.Coordinates(Mesh.Edges(:,2),:))/2;
  Normals = Mesh.Normals;
  
  plot_Mesh(Mesh,'past');
  hold on;
  quiver(MidPts(:,1),MidPts(:,2), ...
         Normals(:,1),Normals(:,2),1/(NREFS+2),'k');
  hold off;
  
  % Clear memory
  
  clear all;
  