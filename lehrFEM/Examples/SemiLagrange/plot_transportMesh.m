function vargout = plot_transportMesh(Mesh,pbVertices,plot_param,varargin)

pMesh.Coordinates=pbVertices(:,[1,2]);
pMesh.Elements=Mesh.Elements;
pMesh = add_Edges(pMesh);
Loc = get_BdEdges(Mesh);
pMesh.BdFlags = zeros(size(pMesh.Edges,1),1);
pMesh.BdFlags(Loc) = -1;

Meshcell{1}=Mesh;
Meshcell{2}=pMesh;

fig=plot_Meshes(Meshcell,plot_param,varargin{:});


% Assign output arguments

if(nargout > 0)
    varargout{1} = fig;
end

return