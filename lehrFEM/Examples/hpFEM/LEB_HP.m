% Refine mesh at each corner nodes and generate HP plot

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    clear
    NRefs = 5;

    % Initialize mesh
    
    Mesh.Coordinates = [0 0;2 0;2 1;3 0;6 3;3 6;0 3;1 2;0 2;3 3];
    Mesh.Elements = [1 8 9;1 3 8;1 2 3;3 4 10;4 5 10;5 6 10;6 7 10;7 8 10;8 3 10];
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags(Loc) = -1;
    Mesh = init_LEB(Mesh);
    CornerNodes = 1:9;
    nCorner = length(CornerNodes);
    
    % Locally refine the mesh
    
    for R = 1:NRefs
        
        Index = [];          % Create container
        for i = 1:nCorner
            [tempI,dummy1,dummy2] = find(Mesh.Elements == CornerNodes(i));
            Index = [Index;tempI];
        end
        Index = unique(Index);
        Mesh = refine_LEB(Mesh,Index);
                
    end
    
    % Assign number to elements

    Mark = tell_Mark(Mesh,CornerNodes,NRefs);
    Mesh.ElemFlag = min(Mark,[],2);
    
    % Plot the mesh with mark
    
    Mesh = rmfield(Mesh,'Edges');
    Mesh = rmfield(Mesh,'BdFlags');
    plot_Mesh(Mesh,'ast');
    
    patch('Faces',Mesh.Elements,'Vertices',Mesh.Coordinates,...
          'FaceVertexCData',Mesh.ElemFlag,'FaceColor','flat');
      
    set(gca,'CLim',[1 max(Mesh.ElemFlag)],'DataAspectRatio',[1 1 1]);
    colormap(jet);
    alpha(.9)
    colorbar
    
    set(gcf,'renderer','openGL')
    print('-depsc', 'LEB_HP.eps');
    !gv LEB_HP.eps &
    
    % Clear memory
    
    clear all
    
    
    