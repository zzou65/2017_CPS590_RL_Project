function arg_out = affine_map(arg_in,Vertices)
% AFFINE_MAP generates the mapping from the reference element 
%
%   COORDINATES = AFFINE_MAP(COORDINATES,VERTICES) generates the points set
%   of each coordinates by the mapping from the reference element to the
%   element which is formed by the given vertices in row-wise orientation
%
%   MESH = AFFINE_MAP(MESH,VERTICES) generates mapping of all the vertices 
%   in the struct MESH by the mapping from the reference element to the
%   element which is formed by the given vertices in row-wise orientation
%
%   VERTICES is a 3-by-2 matrix specifying the vertices of the destination
%   element in a row wise orientation.
%
%   Example:
%
%   Mesh = affine_map(Mesh,[0 0;1 0;1 1]);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Set up element mapping
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
 
  % Map coordinates and assign output arguments
  
  arg_out = deal(arg_in);
  if(isstruct(arg_in))
    nCoordinates = size(arg_in.Coordinates,1);
    arg_out.Coordinates = arg_in.Coordinates*BK+ones(nCoordinates,1)*bK;    
  else
    nCoordinates = size(arg_in,1);
    arg_out = arg_in*BK+ones(nCoordinates,1)*bK;
  end
  
return