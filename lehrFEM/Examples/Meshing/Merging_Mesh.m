% Run script for merging the mesh.

% Copyright 2005-2005 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

    % Initialize constant
    
    NREFS = 1;
    
    % Generate the first mesh

    Mesh_1.Coordinates = [ 0 0 ; 1 0 ; 0 1 ];
    Mesh_1.Elements =[1 2 3];
    Mesh_1 = add_Edges(Mesh_1);
    Loc_1 = get_BdEdges(Mesh_1);
    Mesh_1.BdFlags = zeros(size(Mesh_1.Edges,1),1);
    Mesh_1.BdFlags(Loc_1) = -1;
    for i = 1:NREFS
        Mesh_1=refine_REG(Mesh_1);
    end
    plot_Mesh(Mesh_1,'ptas');

    % Generate the second mesh

    Mesh_2.Coordinates = [ 0 1 ; 1 0 ; 1 1];
    Mesh_2.Elements = [1 2 3];
    Mesh_2 = add_Edges(Mesh_2);
    Loc_2 = get_BdEdges(Mesh_2);
    Mesh_2.BdFlags = zeros(size(Mesh_2.Edges,1),1);
    Mesh_2.BdFlags(Loc_2) = -1;
    for i = 1:NREFS
        Mesh_2=refine_REG(Mesh_2);
    end
    plot_Mesh(Mesh_2,'ptas');

    % Merge the first and second mesh

    Mesh = merge_Mesh(Mesh_1,Mesh_2);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    plot_Mesh(Mesh,'past');

    % Generate the third mesh

    Mesh_3.Coordinates = [ 1 0 ; 2 .5 ; 1 1];
    Mesh_3.Elements = [1 2 3];
    Mesh_3 = add_Edges(Mesh_3);
    Loc_3 = get_BdEdges(Mesh_3);
    Mesh_3.BdFlags = zeros(size(Mesh_3.Edges,1),1);
    Mesh_3.BdFlags(Loc_3) = -1;
    for i = 1:NREFS
        Mesh_3=refine_REG(Mesh_3);
    end
    plot_Mesh(Mesh_3,'past');

    % Merge mesh and mesh3

    Mesh = merge_Mesh(Mesh_3,Mesh);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    plot_Mesh(Mesh,'past');

    % Clear memory

    clear all;
    