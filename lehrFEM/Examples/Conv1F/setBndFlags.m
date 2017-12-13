function Mesh=setBndFlags(Mesh)

Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -2;
for i = Loc'
    if 1/2*sum(Mesh.Coordinates(Mesh.Edges(i,:),2))==0
        Mesh.BdFlags(i) = -1;
    end
end