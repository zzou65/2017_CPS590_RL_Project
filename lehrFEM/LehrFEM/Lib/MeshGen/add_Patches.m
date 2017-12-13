function New_Mesh = add_Patches(Old_Mesh)
% ADD_PATCHES Adds additional patch information to the mesh.
%
%   MESH = ADD_PATCHES(MESH) adds additional patch information to the struct
%   MESH:
%    ADJELEMENTS  M-by-Q matrix specifying the elements of the mesh sharing
%                 vertex i of the mesh.
%    NADJELEMENTS M-by-1 matrix specifying the actual number of elements
%                 sharing vertex COORDINATES(i) of the mesh.  
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the mesh. 
%
%   Example:
%
%   Mesh = add_Patches(Mesh);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   nCoordinates = size(Old_Mesh.Coordinates,1);
   [nElements,nVert] = size(Old_Mesh.Elements);

   % Initialize constants
   
   MAX_ELEMENTS = 10;
   
   % Preallocate memory
   
   AdjElements = zeros(nCoordinates,MAX_ELEMENTS);
   nAdjElements = zeros(nCoordinates,1);
   
   % Build list of adjacent elements
   
   for i = 1:nElements
     for j = 1:nVert
       id = Old_Mesh.Elements(i,j);
       nAdjElements(id) = nAdjElements(id)+1;
       AdjElements(id,nAdjElements(id)) = i;
     end
   end
   
   % Squeeze arrays
   
   MAX_ELEMENTS = max(nAdjElements);
   AdjElements = AdjElements(:,1:MAX_ELEMENTS);
   
   % Assign output arguments
   
   New_Mesh = deal(Old_Mesh);
   New_Mesh.AdjElements = AdjElements;
   New_Mesh.nAdjElements = nAdjElements;
   
return
