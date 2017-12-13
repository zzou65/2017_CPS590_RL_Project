function Mesh = add_Elements(Mesh)
% ADD_ELEMENTS Add elements list
%
%   MESH = ADD_ELEMENTS(MESH) adds elements to a mesh using the delaunay 
%   algorithm
%
%   The struct MESH must contain:
%    COORDINATES M-by-2 matirx specifying the vertices of the mesh.
%
%   Example:
%
%   x = linspace(-1,1,10);
%   y = linspace(-1,1,10);
%   [X Y] = meshgrid(x,y);
%   Mesh.Coordinates = [X(:,1) X(:,2)];
%   Mesh = add_Elements(Mesh);

%   Copyright 2005-2005 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Add elements

  try
      Mesh.Elements = delaunayn(Mesh.Coordinates);
  catch
      Mesh.Elements = delaunay(Mesh.Coordinates(:,1),Mesh.Coordinates(:,2),{'Qt','Qbb','Qc','Qz'});
  end

  % Change orientation of all elements to counter-clockwise
  
  for i = 1:size(Mesh.Elements,1);
    v = Mesh.Coordinates(Mesh.Elements(i,2),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    w = Mesh.Coordinates(Mesh.Elements(i,3),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    if(v(1)*w(2)-v(2)*w(1) < 0)
      tmp = Mesh.Elements(i,1);
      Mesh.Elements(i,1) = Mesh.Elements(i,2);
      Mesh.Elements(i,2) = tmp;
    end
  end
  
  return