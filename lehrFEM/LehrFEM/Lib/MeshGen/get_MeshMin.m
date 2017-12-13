function h = get_MeshMin(Mesh)
% GET_MESHMIN Computes the smallest edge length
%   
%   H = MESH_MESHMIN(MESH) computes the minimal length of edges for a given
%   triangular or quarilateral mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of
%                the mesh.
%
%   Example:
%
%   h = get_MeshMin(Mesh);

%   Copyright 2005-2006 Patrick Meury and Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   % Intialize constants

   nVert = size(Mesh.Elements,2);
 
   if(nVert == 3)
   
    % Compute mesh width for triangular meshes   
       
    h1 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2)));
    h2 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2)));
    h3 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,2),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,1),:)).^2,2)));
    h = min([h1 h2 h3]);
    
   else
       
     % Compute mesh width for  quadrilateral meshes  
       
     h1 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,2),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,1),:)).^2,2)));
     h2 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2)));
     h3 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,4),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2)));
     h4 = min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:) - ...
                       Mesh.Coordinates(Mesh.Elements(:,4),:)).^2,2)));
     h = min([h1 h2 h3 h4]);
     
   end
    
return
    
           