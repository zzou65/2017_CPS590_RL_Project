function Mesh=create_Mesh(N, sigma)

% Mesh_1=createMesh_1(0,0,2);
% plot_Mesh(Mesh_1,'petas')
% 
% Mesh_2=createMesh_2(0,0,1);
% plot_Mesh(Mesh_2,'petas')
% 
% Mesh_3=createMesh_3(0,0,1);
% plot_Mesh(Mesh_3,'petas')
% 
% Mesh_4=createMesh_4(0,0,1);
% plot_Mesh(Mesh_4,'petas')
% 
%  Example Mesh=create_Mesh(10)

npoints = 2.^floor(N*sigma);
points=[1/npoints:1/npoints:1-1/npoints];
N=2^N;
h=1/N;
Mesh=createMesh_1(h/4,h/4,h);
for i=1:(N-1)
    Mesh_2=createMesh_2((i-1/2)*h,h/4,h);
    % merge
    Mesh = merge_Mesh(Mesh_2,Mesh);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    
    Mesh_3=createMesh_3((i)*h,h/4,h);
    % merge
    Mesh = merge_Mesh(Mesh_3,Mesh);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
end
Mesh_2=createMesh_2((N-1/2)*h,h/4,h);
% merge
Mesh = merge_Mesh(Mesh_2,Mesh);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

Mesh_4=createMesh_4((N-1/4)*h,h/4,h);
% merge
Mesh = merge_Mesh(Mesh_4,Mesh);
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;

FlippedMesh = Mesh;
FlippedMesh.Coordinates(:,2) = h-FlippedMesh.Coordinates(:,2);
FlippedMesh = add_Edges(FlippedMesh);
Loc = get_BdEdges(FlippedMesh);
FlippedMesh.BdFlags = zeros(size(FlippedMesh.Edges,1),1);
FlippedMesh.BdFlags(Loc) = -1;
FlippedMesh = orient_Elems(FlippedMesh);
% merge
MeshLine = merge_Mesh(FlippedMesh,Mesh);
MeshLine = add_Edges(MeshLine);
Loc = get_BdEdges(MeshLine);
MeshLine.BdFlags = zeros(size(MeshLine.Edges,1),1);
MeshLine.BdFlags(Loc) = -1;

Mesh = MeshLine;
for i=1:(N-1)
    MeshLine.Coordinates(:,2) = MeshLine.Coordinates(:,2)+h;
    % merge
    Mesh = merge_Mesh(MeshLine,Mesh);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
end

Mesh = init_LEB(Mesh);
if ~isempty(points)
    Melements=markElements(Mesh,[points]);
    Mesh = refine_LEB(Mesh,Melements);
    % Update mesh data structure
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;
end
%plot_Mesh(Mesh,'as')


function Mesh_1=createMesh_1(x0,y0,h)
Mesh_1.Coordinates = [ x0-0.25*h y0-0.25*h ; x0+0.25*h y0+0.25*h ; x0-0.25*h y0+0.25*h];
Mesh_1.Elements =[1 2 3];
Mesh_1 = add_Edges(Mesh_1);
Loc_1 = get_BdEdges(Mesh_1);
Mesh_1.BdFlags = zeros(size(Mesh_1.Edges,1),1);
Mesh_1.BdFlags(Loc_1) = -1;

function Mesh_2=createMesh_2(x0,y0,h)
Mesh_2.Coordinates = [ x0+0 y0+0.25*h ; x0-0.5*h y0-0.25*h ; x0+0.5*h y0-0.25*h];
Mesh_2.Elements =[1 2 3];
Mesh_2 = add_Edges(Mesh_2);
Loc_2 = get_BdEdges(Mesh_2);
Mesh_2.BdFlags = zeros(size(Mesh_2.Edges,1),1);
Mesh_2.BdFlags(Loc_2) = -1;

function Mesh_3=createMesh_3(x0,y0,h)
Mesh_3.Coordinates = [ x0+0 y0-0.25*h ; x0+0.5*h y0+0.25*h ; x0-0.5*h y0+0.25*h];
Mesh_3.Elements =[1 2 3];
Mesh_3 = add_Edges(Mesh_3);
Loc_3 = get_BdEdges(Mesh_3);
Mesh_3.BdFlags = zeros(size(Mesh_3.Edges,1),1);
Mesh_3.BdFlags(Loc_3) = -1;

function Mesh_4=createMesh_4(x0,y0,h)
Mesh_4.Coordinates = [ x0-0.25*h y0+0.25*h ; x0+0.25*h y0-0.25*h ; x0+0.25*h y0+0.25*h];
Mesh_4.Elements =[1 2 3];
Mesh_4 = add_Edges(Mesh_4);
Loc_4 = get_BdEdges(Mesh_4);
Mesh_4.BdFlags = zeros(size(Mesh_4.Edges,1),1);
Mesh_4.BdFlags(Loc_4) = -1;

function MarkedElements=markElements(Mesh,ypoints)

m=length(ypoints);
nElements= size(Mesh.Elements,1);
MarkedElements=[];

for i=1:nElements
    if ~isempty(intersect(ypoints,1/3*sum(Mesh.Coordinates(Mesh.Elements(i,:),1))));
        MarkedElements=[MarkedElements,i];
    end
end

