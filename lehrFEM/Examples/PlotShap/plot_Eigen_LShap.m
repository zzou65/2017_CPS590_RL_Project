    
    % Clear memory
    
    clear all;
    
    % Initialize constant

    NEigen = 3;                          % Number of the eigenvalue
    NREFS = 5;                            % Number of red refinement steps

    % Initialize mesh
    
    Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
    LShap_Vert = Mesh.Coordinates;              % Extract outline
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    for i = 1:NREFS
        Mesh = refine_REG(Mesh);
    end   
    
    % Assemble stiffness matrix and mass matrix
    
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
    M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    
    % Solve eigenvalue problem
    
    [ V d ] = eigs(A,M,NEigen,'sm');
    
    for i = 1:NEigen
        
        U(:,i) = V(:,i);
        norm(A*V(:,i)-d(i,i)*M*V(:,i))
      
    end

    % Plot eigen function in 3D

    % Initialize constant

    p = 1;                 % Background color

    % Generate mesh

    for i = 1:size(Mesh.Elements,1);
        v = Mesh.Coordinates(Mesh.Elements(i,2),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
        w = Mesh.Coordinates(Mesh.Elements(i,3),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
        if(v(1)*w(2)-v(2)*w(1) > 0)
            tmp = Mesh.Elements(i,1);
            Mesh.Elements(i,1) = Mesh.Elements(i,2);
            Mesh.Elements(i,2) = tmp;
        end
    end

    % Adjust solution
    
    U(:,2) = -U(:,2);
    U = U(:,3:-1:1);

    % Plot shap function

    OFFSET = 0.05;
    
    for i=2:NEigen

        % Output outline
        
        plot3([LShap_Vert(:,1);LShap_Vert(1,1)], ...
              [LShap_Vert(:,2);LShap_Vert(1,2)], ...
              zeros(length(LShap_Vert)+1,1),'Color',[.5 .5 .5],'linewidth',1.5)
          
        % Compute axes limits

        XMin = min(Mesh.Coordinates(:,1));
        XMax = max(Mesh.Coordinates(:,1));
        YMin = min(Mesh.Coordinates(:,2));
        YMax = max(Mesh.Coordinates(:,2));
        XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
        YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];

        % Compute color axes limits

        CMin = min(U(:,i));
        CMax = max(U(:,i));
        if(CMin < CMax)          % or error will occur in set function
            CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
        else
            CLim = [1-OFFSET 1+OFFSET]*CMin;
        end
        set(gcf,'Color',[1 1 1]*p);
        set(gca,'CameraPosition', [1 -.5 .5]*50,...
            'CameraTarget',[0 0 0], ...
            'CameraUpVector',[0 0 1], ...
            'CameraViewAngle',5, ...
            'DataAspectRatio', [1 1 0.9],...
            'Position',[0 .1 1 1], ...
            'LineWidth',2,...
            'Visible','on', ...
            'XLim',[-1.1 1.1], ...
            'YLim',[-1.1 1.1], ...
            'ZLim',[CMin,CMax]);
        xlabel('\bfX')
        ylabel('\bfY')
        zlabel('\bfZ')
        patch('faces', Mesh.Elements, ...
            'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U(:,i)], ...
            'CData', U, ...
            'EdgeColor','none', ...
            'FaceColor',[0.9 0.2 0.2], ...
            'FaceLighting','phong', ...
            'AmbientStrength',0.3, ...
            'DiffuseStrength',0.6, ...
            'Clipping','off',...
            'BackFaceLighting','lit', ...
            'SpecularStrength',.9, ...
            'SpecularColorReflectance',1, ...
            'SpecularExponent',7);
        alpha(.7)

        l1 = light('Position',[-1 -.5 1.3], ...
            'Style','infinite', ...
            'Color',[0 0.8 0.8]);
        l2 = light('Position',[1 1 2], ...
            'Style','local', ...
            'Color',[0.8 .8 0]);
        l3 = light('Position',[1 -1 .5], ...
            'Style','local', ...
            'Color',[0.8 .8 0]);

        FileName = ['Eigen_LShap3D_' int2str(i-1) '.eps'];
        print('-depsc', FileName);
        system(['gv ' FileName ' &']);

    end
    
    % Clear memory
    
    clear all;
    
    