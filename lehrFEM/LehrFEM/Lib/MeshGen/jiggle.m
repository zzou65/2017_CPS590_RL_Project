function Mesh = jiggle(Mesh,FixedPos)
% JIGGLE Mesh Jiggling.
%
%   MESH = JIGGLE(MESH,FIXEDPOS) jiggles the struct MESH using the 
%   jiggling algorithm.
%
%   The array FIXEDPOS specifies the fixed vertices of the mesh:
%    0 Vertex position will be moved.
%    1 Vertex position will not be moved.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements
%                of the mesh. 
%
%   Example:
%
%   Mesh = jiggle(Mesh,FixedPos);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  FixedPos = find(FixedPos);

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);  % Number of vertices
  nFixedPos = size(FixedPos,1);             % Number of fixed vertices
  nIntPos = nCoordinates-nFixedPos;         % Number of interior vertices
  sca = 0.2;                                % Shift scaling parameter           
  
  if(size(Mesh.Elements,2) == 3)
  
    % Compute minimal edge length for triangles
  
    h = min([min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2))) ...
             min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2))) ...
             min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2)))]);
                         
  else
  
    % Compute minimal edge length quadrilaterals
  
    h = min([min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,1),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,2),:)).^2,2))) ...
             min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,2),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,3),:)).^2,2))) ...
             min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,3),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,4),:)).^2,2))) ...
             min(sqrt(sum((Mesh.Coordinates(Mesh.Elements(:,4),:)- ...
                           Mesh.Coordinates(Mesh.Elements(:,1),:)).^2,2)))]);    
   
  end                
                   
  % Initialize random shift vectors
   
  theta = 2*pi*rand(nIntPos,1);
  r = sca*h*rand(nIntPos,1);
  shift = [r.*cos(theta) r.*sin(theta)];
   
  % Shift all interior vertices of the mesh
     
  IntPos = setdiff(1:nCoordinates,FixedPos);
  Mesh.Coordinates(IntPos,:) = Mesh.Coordinates(IntPos,:) + shift;  
      
return
   