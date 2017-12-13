% Run script for generating graded mesh of L-Shape domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


    % Initialize constant
    
    Beta = 4;                             % Grading factor
    nNodes = 65;                          % Number of nodes
    FHandle = @(x,varagin)x.^Beta;        % Algebraically graded
    FHandle2 = @(x,varagin)x;             % Algebraically graded 2
    
    % Generate each part of L-Shape by mesh mapping
    
    % Generate graded mesh for the reference element
    
    Mesh_REF = graded_RefElem(nNodes,FHandle);

    Mesh_1 = affine_map(Mesh_REF,[ 0 0 ; 1 0 ; 0 1]); 
    Mesh_1.Elements = delaunayn(Mesh_1.Coordinates);
    Mesh_1 = orient_Elems(Mesh_1);
    Mesh_1 = add_Edges(Mesh_1);
    
    Mesh_2 = affine_map(Mesh_REF,[ 0 0 ; 0 1 ; -1 0]);
    Mesh_2.Elements = delaunayn(Mesh_2.Coordinates);
    Mesh_2 = orient_Elems(Mesh_2);
    Mesh_2 = add_Edges(Mesh_2);
  
    Mesh_3 = affine_map(Mesh_REF,[ 0 0 ; -1 0 ; 0 -1]);
    Mesh_3.Elements = delaunayn(Mesh_3.Coordinates);
    Mesh_3 = orient_Elems(Mesh_3);  
    Mesh_3 = add_Edges(Mesh_3);

    % Generate part2
    
    Mesh_REF2 = graded_RefElem(nNodes,FHandle2);
    
    Mesh_4 = affine_map(Mesh_REF2,[ 1 1 ; 0 1 ; 1 0]);
    Mesh_4.Elements = delaunayn(Mesh_4.Coordinates);
    Mesh_4 = orient_Elems(Mesh_4);
    Mesh_4 = add_Edges(Mesh_4);
    
    Mesh_5 = affine_map(Mesh_REF2,[ -1 1 ; -1 0 ; 0 1]);
    Mesh_5.Elements = delaunayn(Mesh_5.Coordinates);
    Mesh_5 = orient_Elems(Mesh_5);
    Mesh_5 = add_Edges(Mesh_5);
    
    Mesh_6 = affine_map(Mesh_REF2,[ -1 -1 ; 0 -1; -1 0]);
    Mesh_6.Elements = delaunayn(Mesh_6.Coordinates);
    Mesh_6 = orient_Elems(Mesh_6);
    Mesh_6 = add_Edges(Mesh_6);
    
    % Merge the meshes
    
    Mesh = merge_Mesh(Mesh_1,Mesh_2);
    Mesh = merge_Mesh(Mesh,Mesh_3);    
    Mesh = merge_Mesh(Mesh,Mesh_4);
    Mesh = merge_Mesh(Mesh,Mesh_5);
    Mesh = merge_Mesh(Mesh,Mesh_6); 
     
    % Plot mesh
    
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    plot_Mesh(Mesh,'as');
    plot_Qual(Mesh);
    
    
    % Clear Memory
    
    clear all;
     