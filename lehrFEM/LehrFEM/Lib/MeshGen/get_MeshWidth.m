function h = get_MeshWidth(Mesh)
% GET_MESHWIDTH Computes the mesh width.
%   
%   H = MESH_MESHWIDTH(MESH) computes the mesh width for a given triangular or
%   quarilateral mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of
%                the mesh.
%
%   Example:
%
%   h = get_MeshWidth(Mesh);

%   Copyright 2005-2005 Patrick Meury and Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   % Intialize constants

   nVert = size(Mesh.Elements,2);
 
   if(nVert == 3)
   
    % Compute mesh width for triangular meshes   
       
    h1 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2)));
    h2 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2)));
    h3 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,2),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,1),:)).^2,2)));
    h = max([h1 h2 h3]);
    
   else
       
     % Compute mesh width for  quadrilateral meshes  
       
     h1 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,2),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,1),:)).^2,2)));
     h2 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2)));
     h3 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,4),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2)));
     h4 = max(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,4),:)).^2,2)));
     h = max([h1 h2 h3 h4]);
     
   end
    
return
    
                          