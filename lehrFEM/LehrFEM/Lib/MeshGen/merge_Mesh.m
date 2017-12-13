function New_Mesh = merge_Mesh(Mesh_1,Mesh_2)
% MERGE_MESH merge Mesh2 to Mesh1 and generate a new mesh
%
%   MESH = MERGE_MESH(MESH_1,MESH_2) merge MESH_1 to MESH_2 and generate a 
%   new mesh.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   Notice: MESH_1 and MESH_2 must be matched on the interface
%
%   Example:
%
%   Mesh = merge_Mesh(Mesh_1,Mesh_2);
%
%   See also add_Edges, get_BdEdges.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang 
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

   % Exchange meshes if the second one is larger than the first one
   
   if (size(Mesh_1.Coordinates,1)<size(Mesh_2.Coordinates,1))
       TempMesh = Mesh_2;
       Mesh_2 = Mesh_1;
       Mesh_1 = TempMesh;
   end
   
   % Initialize constants
   
   TOL = eps;
   nCoordinates_1 = size(Mesh_1.Coordinates,1);
   nElements_1 = size(Mesh_1.Elements,1);
   nCoordinates_2 = size(Mesh_2.Coordinates,1);
   nElements_2 = size(Mesh_2.Elements,1);
   
   % Precompute mesh information
   
   Loc_1 = get_BdEdges(Mesh_1);
   BdNodes_1 = unique([Mesh_1.Edges(Loc_1,1) Mesh_1.Edges(Loc_1,2)]); 
   nBdNodes_1 = length(BdNodes_1);
   V1 = Mesh_1.Coordinates(BdNodes_1,:);
   
   Loc_2 = get_BdEdges(Mesh_2);
   BdNodes_2 = unique([Mesh_2.Edges(Loc_2,1) Mesh_2.Edges(Loc_2,2)]); 
   nBdNodes_2 = length(BdNodes_2);
   V2 = Mesh_2.Coordinates(BdNodes_2,:);
   
   % Preallocate memory 
   
   Match = zeros(nBdNodes_1,1);
   VertMap = zeros(nCoordinates_2,1);
   
   % Find joint points between two meshes
   
   for i = 1:nBdNodes_1
       res = V2 - ones(nBdNodes_2,1)*V1(i,:);
       res = sqrt(sum(res.^2,2));
       Mark = find(abs(res)<TOL);
       if(~isempty(Mark))
           Match(i) = Mark;
       end
   end
   
   % Built vertex mapping
   
   [Joint_1 dummy Joint_2] = find(Match);
   nJointNodes = length(Joint_2);
   Joint_1 = BdNodes_1(Joint_1);
   Joint_2 = BdNodes_2(Joint_2);  
   VertMap(Joint_2) = Joint_1;
   Non_Joint2 = setdiff(1:nCoordinates_2,Joint_2);
   VertMap(Non_Joint2) = (1: nCoordinates_2-nJointNodes) + nCoordinates_1;
   
   % Rearrange mesh2
   
   Mesh_2.Elements(:,1) = VertMap(Mesh_2.Elements(:,1));
   Mesh_2.Elements(:,2) = VertMap(Mesh_2.Elements(:,2));
   Mesh_2.Elements(:,3) = VertMap(Mesh_2.Elements(:,3));
   
   % Generate new mesh
   
    % Initialize constant
    
    nCoordinates  = nCoordinates_1+nCoordinates_2-nJointNodes;
    nElements = nElements_1+nElements_2;
   
    % Preallocate memory
    
    New_Mesh.Coordinates = zeros(nCoordinates,2);
    New_Mesh.Elements = zeros(nElements,3);
    
    % Generate coordinates field to the new mesh
       
    New_Mesh.Coordinates = [ Mesh_1.Coordinates ; Mesh_2.Coordinates(Non_Joint2,:) ];
   
    % Generate Elements and Edges fields to the new mesh
   
    New_Mesh.Elements = [ Mesh_1.Elements ; Mesh_2.Elements];
   
    New_Mesh = add_Edges(New_Mesh);
    
return
  