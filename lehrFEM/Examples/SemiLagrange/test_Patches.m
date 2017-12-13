% test the function patches and determine the elements where the sum of the
% area of the patches does not match the total element size
% do a plot with the errors for each element
function test2(NewMesh,defMesh,intersec); 


nElements=size(NewMesh.Elements,1);
err = zeros(nElements,1);

% test the patches
totarea = 0;
for i=1:nElements
    nPatches = intersec(i).nElems;
    polygons = patches(NewMesh,defMesh,i,intersec(i));
    %i
    area = 0;
    vert = NewMesh.Elements(i,:);
    e1 = NewMesh.Coordinates(vert(1),:);
    e2 = NewMesh.Coordinates(vert(2),:);
    e3 = NewMesh.Coordinates(vert(3),:);
    m = [e1;e2;e3];
    elemarea = polyarea(m(:,1),m(:,2));
    for j=1:nPatches
        polygon = polygons(j).Polygon;
        if ~isempty(polygon)
            area = area + polyarea(polygon(:,1),polygon(:,2));
        end
    end
    
    totarea = totarea + area;
    err(i) = abs(elemarea - area);
    
    if abs(elemarea - area) > 10^-4
        i
        elemarea - area
    end
end
plot_P0(err,NewMesh); colorbar;

totarea

%addpath([pwd '/Examples/SemiLagrange/']);