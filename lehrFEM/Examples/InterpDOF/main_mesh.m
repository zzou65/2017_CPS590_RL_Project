% Run script for interest

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constant
    
    NREFS =7;
    lambda = 0.2:0.2:2.8;
    %lambda = linspace(1.9,2.1,5);
    lambda = [0.2:0.2:1.8 1.9 2.2:0.2:2.8];
    NLAMBDA = length(lambda);
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,NLAMBDA);
    Err_H1 = zeros(NREFS,NLAMBDA);
    p = zeros(NLAMBDA,2);
    % Initialize mesh

    X0 = [-1 -1];                                     % Lower left corner point of rectangle
    A = 2;                                          % Length of rectangle
    BBOX = [X0; X0+[A A]];                          % Bounding box
    H0 = 0.2;                                       % Initial mesh width
    DHANDLE = @dist_rect;                           % Signed distance function
    HHANDLE = @h_uniform;                           % Element size function
    FIXEDPOS = [X0; X0+[A 0]; X0+[A A]; X0+[0 A]];  % Fixed boundary vertices of the mesh
    DISP = 0;                                       % Display flag
    Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,X0,A,A);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    
    % Run the script
    
    for iter = 1:NREFS
        
        % Do REG refine
        Mesh = refine_REG_jiggle(Mesh);
        %Mesh = refine_REG(Mesh);
        figure;
        plot_Mesh(Mesh,'as');
        file = ['Mesh_' int2str(iter) '.eps'];
        print('-depsc',file);
    end
    
    clear all;