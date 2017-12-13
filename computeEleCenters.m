function eleCenterCoords = computeEleCenters(Mesh)
 
nElements = size(Mesh.Elements,1);
eleCenterCoords = zeros(nElements, 2);
for i=1:nElements
    node_idx = Mesh.Elements(i,:);
    Vertices_Coors = Mesh.Coordinates(node_idx,:);
    eleCenterCoords(i, :) = mean(Vertices_Coors);
end

end