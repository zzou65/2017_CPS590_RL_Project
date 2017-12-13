function [mesh] = Sqr_Mesh(R)

mesh = [];
mesh.Coordinates = [0 0;1 0;1 1;0 1];
mesh.Elements = [1 2 3;1 3 4];
mesh = add_Edges(mesh);
mesh.BdFlags = [-1;0;-1;-1;-1];

for i = 1:R
    mesh = refine_REG(mesh);
end;
