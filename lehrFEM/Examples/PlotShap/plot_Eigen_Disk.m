% Run script for the eigenvalue problem in the unit square domain

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    NEigen = 2;                    % Number of the eigenvalue
    Select = [1 2 4 6 7 9 11 13 15 16 18 20 22 24];
    ZBessel = [2.4048 3.8317 5.1356 6.3802 7.5883 8.7715 5.5201...
              7.0156 8.4172 9.7610 11.0647 8.6537 10.1735 9.93611];
    ZBessel = sort(ZBessel)';
    L = 1:14;
    Lambda = ZBessel.^2;
    
    % Initialize mesh
    
    C = [0 0];             % Center of circle
    R = 1;                 % Radius of circle 
    BBOX = [-1 -1; 1 1];   % Bounding box
    H0 = 0.05;             % Initial mesh width
    DHANDLE = @dist_circ;  % Signed distance function
    HHANDLE = @h_uniform;  % Element size function
    FIXEDPOS = [];         % Fixed boundary vertices of the mesh
    DISP = 0;              % Display flag

    Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);  
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
 
    % Assemble stiffness matrix and mass matrix
    
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
    M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    
    
    % Incorporate Neumann boundary data
    
    [U,FreeNodes] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    A = A(FreeNodes,FreeNodes);
    M = M(FreeNodes,FreeNodes);

    % Solve eigenvalue problem
    
    [ V d ] = eigs(A,M,NEigen,'sm');
    
    for i = 1:NEigen
        
        U(FreeNodes,i) = V(:,i);
        norm(A*V(:,i)-d(i,i)*M*V(:,i))
      
    end

    % Plot eigen function in 3D

    % Initialize constant

    p = .9;                 % Background color
    
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

    % Generate outline
    
    t = linspace(0,2*pi,200);
    tx = sin(t);
    ty = cos(t);
    
    % Adjust solution
    
    U(:,2) = -U(:,2);

    % Plot shap function

    OFFSET = 0.05;

    for i=1:NEigen

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

        % Plot outline
        
        fig = figure;
        plot3(tx,ty,zeros(size(tx)),'-k','linewidth',1.5)
        
        % Generate figure
        
        set(gcf,'Color',[1 1 1]*p);
        set(gca,'CameraPosition', [-1 -1 .9]*50,...
            'CameraTarget',[.25 .25 0], ...
            'CameraUpVector',[0 0 1], ...
            'CameraViewAngle',3.5, ...
            'DataAspectRatio', [1 1 0.9],...
            'Position',[0 0.1 1 1], ...
            'Visible','on', ...
            'LineWidth',2,...
            'XLim',XLim, ...
            'YLim',YLim, ...
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
        
        l1 = light('Position',[-1 -.5 -.3], ...
            'Style','infinite', ...
            'Color',[0 0.8 0.8]);
        l2 = light('Position',[7 7 10], ...
            'Style','infinite', ...
            'Color',[0.8 .8 0]);

        FileName = ['Eigen_Disk3D_' int2str(NEigen-i+1) '.eps'];
        print('-depsc', FileName);
        system(['gv ' FileName ' &']);

    end