
% Create Mesh
Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
Mesh = add_Edges(Mesh);
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
Mesh = refine_REG(Mesh);
Mesh = refine_REG(Mesh);
Mesh=add_Edge2Elem(Mesh);

% Velocityfield
v1=0.4;
v2=0.2;
Dir_Handle=@(x,varargin) -(1-x(:,2).^2).*(1-x(:,1).^2)*[v1 v2];

% transport of Coordinates
directions = Dir_Handle(Mesh.Coordinates);
dt=0.1*get_MeshWidth(Mesh);
pulled_back_vertices=trace_vertices(Mesh,-dt*directions);

% plot transported Mesh
plot_transportMesh(Mesh,pulled_back_vertices,'petas');

% stiffness matrices 

% assemble by coordinates, precompute weights
M_0 =  assemMat_Mass0fD(Mesh);
A =assemMat_SemiLag_LFE(Mesh,pulled_back_vertices);
% assemble by elements
B =assemMat_SemiLagQuad_LFE(Mesh,pulled_back_vertices); 

% It should be M_0 * A == B;
A
B


