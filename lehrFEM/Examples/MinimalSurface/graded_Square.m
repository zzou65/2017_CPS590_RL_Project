function Mesh = graded_Square(nNodes,Beta,varargin)

% GRADED_SQUARE generates a graded mesh of the unit square
%
%   GRADED_SQUARE(NNODES,BETA) generates a graded square mesh of the square
%   with corner coordinates (-1,-1),(1,-1),(1,1),(-1,1) with grading factor
%   BETA. The final mesh is merged together by 8 triangular meshes, each 
%   having NNODES unknowns along each katet. The mesh is graded towards the
%   points (-1,0) and (1,0).
%
%   NNODES is a measurement of the total number of unknowns
%   
%   BETA is the grading factor
%
%   GRADED_SQUARE(NNODES,BETA,'P') also plots the mesh and the quality
%   plot of the mesh
%   
%   Example
%
%   Mesh = graded_Square(10,2);

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland



    FHandle = @(x,varagin)x.^Beta;  
    FHandle2 = @(x,varagin)x;            
        
    % Generate graded mesh for the reference element
    
    Mesh_REF = graded_RefElem(nNodes,FHandle);

    Mesh_1 = affine_map(Mesh_REF,[-1 0 ; 0 0 ; -1 1]); 
    Mesh_1.Elements = delaunayn(Mesh_1.Coordinates);
    Mesh_1 = orient_Elems(Mesh_1);
    Mesh_1 = add_Edges(Mesh_1);
    
    Mesh_2 = affine_map(Mesh_REF,[-1 0 ; -1 -1 ; 0 0]);
    Mesh_2.Elements = delaunayn(Mesh_2.Coordinates);
    Mesh_2 = orient_Elems(Mesh_2);
    Mesh_2 = add_Edges(Mesh_2);
  
    Mesh_3 = affine_map(Mesh_REF,[ 1 0 ; 1 1 ; 0 0]);
    Mesh_3.Elements = delaunayn(Mesh_3.Coordinates);
    Mesh_3 = orient_Elems(Mesh_3);  
    Mesh_3 = add_Edges(Mesh_3);
    
    Mesh_4 = affine_map(Mesh_REF,[ 1 0 ; 0 0 ; 1 -1]);
    Mesh_4.Elements = delaunayn(Mesh_4.Coordinates);
    Mesh_4 = orient_Elems(Mesh_4);  
    Mesh_4 = add_Edges(Mesh_4);

    % Generate part2
    
    Mesh_REF2 = graded_RefElem(nNodes,FHandle2);
    
    Mesh_5 = affine_map(Mesh_REF2,[0 0 ; 0 1 ; -1 1]);
    Mesh_5.Elements = delaunayn(Mesh_5.Coordinates);
    Mesh_5 = orient_Elems(Mesh_5);
    Mesh_5 = add_Edges(Mesh_5);
    
    Mesh_6 = affine_map(Mesh_REF2,[0 0 ; 1 1 ; 0 1]);
    Mesh_6.Elements = delaunayn(Mesh_6.Coordinates);
    Mesh_6 = orient_Elems(Mesh_6);
    Mesh_6 = add_Edges(Mesh_6);
    
    Mesh_7 = affine_map(Mesh_REF2,[0 0 ; -1 -1; 0 -1]);
    Mesh_7.Elements = delaunayn(Mesh_7.Coordinates);
    Mesh_7 = orient_Elems(Mesh_7);
    Mesh_7 = add_Edges(Mesh_7);
    
    Mesh_8 = affine_map(Mesh_REF2,[0 0 ; 0 -1; 1 -1]);
    Mesh_8.Elements = delaunayn(Mesh_8.Coordinates);
    Mesh_8 = orient_Elems(Mesh_8);
    Mesh_8 = add_Edges(Mesh_8);
    
    % Merge the meshes
    
    Mesh = merge_Mesh(Mesh_1,Mesh_2);
    Mesh = merge_Mesh(Mesh,Mesh_3);    
    Mesh = merge_Mesh(Mesh,Mesh_4);
    Mesh = merge_Mesh(Mesh,Mesh_5);
    Mesh = merge_Mesh(Mesh,Mesh_6);
    Mesh = merge_Mesh(Mesh,Mesh_7);
    Mesh = merge_Mesh(Mesh,Mesh_8);
    
    if nargin == 3,
        
        % Plot mesh
        
        Loc = get_BdEdges(Mesh);
        Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
        Mesh.BdFlags(Loc) = -1;
        plot_Mesh(Mesh,'as');
        plot_Qual(Mesh);
    end
        