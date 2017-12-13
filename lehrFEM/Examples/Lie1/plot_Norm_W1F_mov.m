function PU = plot_Norm_W1_movF(U,Mesh)
% PLOT_NORM_W1F Plot routine for the modulus of W1F results.
%
%   FIG = PLOT_NORM_W1F(U,MESH) generates a plot of the modulus for the 
%   velocity field which is represented by the W1F solution U on the struct
%   MESH and returns the handle FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh.
%    EDGES       P-by-2 matrix specifying the edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies whether the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%
%   Example:
%
%   fig = plot_Norm_W1F(U,Mesh);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    nElements = size(Mesh.Elements,1);
    nCoordinates = size(Mesh.Coordinates,1);
    
    % Preallocate memory
    
    ux = zeros(nElements,1);
    uy = zeros(nElements,1);
    PU = zeros(nCoordinates,1);
    PUx = zeros(nCoordinates,1);
    PUy = zeros(nCoordinates,1);
    
    % Calculate modulus
    
    for i = 1:nElements
        
        vidx = Mesh.Elements(i,:);
        P1 = Mesh.Coordinates(vidx(1),:);
        P2 = Mesh.Coordinates(vidx(2),:);
        P3 = Mesh.Coordinates(vidx(3),:);
        
        bK = P1;
        BK = [P2-P1;P3-P1];
        TK = transpose(inv(BK));
       
        % Locate barycenter
        
        Bar_Node = 1/3*[P1 + P2 + P3];
        
        % Compute velocity field at barycenters

        eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
                Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
                Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];
        
        % Determine edge orientation

        if(Mesh.Edges(eidx(1),1) == vidx(2))
            p1 = 1;
        else
            p1 = -1;
        end

        if(Mesh.Edges(eidx(2),1) == vidx(3))
            p2 = 1;
        else
            p2 = -1;
        end

        if(Mesh.Edges(eidx(3),1) == vidx(1))
            p3 = 1;
        else
            p3 = -1;
        end
        
        % Compute velocity field at barycenters
        
        N = shap_W1F(Bar_Node);
        NS(1:2) = N(1:2)*TK;
        NS(3:4) = N(3:4)*TK;
        NS(5:6) = N(5:6)*TK;
        
        ux(i) = U(eidx(1))*p1*NS(1) + ...
                U(eidx(2))*p2*NS(3) + ...
                U(eidx(3))*p3*NS(5);
        
        uy(i) = U(eidx(1))*p1*NS(2) + ...
                U(eidx(2))*p2*NS(4) + ...
                U(eidx(3))*p3*NS(6);
                 
    end
    
    % Calculate value on each vertice
    
    Mesh = add_Patches(Mesh);
    
    for i = 1:nCoordinates
        L_patch = Mesh.AdjElements(i,:);
        loc = find(L_patch>0);
        Eidx = L_patch(loc);
        PUx(i) = sum(ux(Eidx))/Mesh.nAdjElements(i);   
        PUy(i) = sum(uy(Eidx))/Mesh.nAdjElements(i);   
    end

    PU = sqrt(PUx.^2+PUy.^2);
    
return

