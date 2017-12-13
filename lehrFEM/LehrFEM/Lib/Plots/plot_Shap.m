function varargout = plot_Shap(vertex,U)
% PLOT_SHAP plot U on the set of vertex in 3D with lighting effect
%
%   PLOT_SHAP(VERTEX,U) generates a 3D plot of U on the domain represented 
%   by vertex
%   
%   The matrix Vertex determine the set of point which represents the
%   domain of plotting.
%
%   H = PLOT_SHAP(VERTEX,U) also returns the handle to the 
%   figure.
%
%   Example:
%
%   plot_Shap([0 0;1 0;0 1],[1;2;3]);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


    % Initialize constant
    
    p = 1;                  % Background color
    
    % Generate mesh
    
    elem = delaunayn(vertex);
    Mesh.Coordinates = vertex;
    Mesh.Elements = elem;
    
    % Rearrange orientations
    
    for i = 1:size(Mesh.Elements,1);
        v = Mesh.Coordinates(Mesh.Elements(i,2),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
        w = Mesh.Coordinates(Mesh.Elements(i,3),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
        if(v(1)*w(2)-v(2)*w(1) > 0)
            tmp = Mesh.Elements(i,1);
            Mesh.Elements(i,1) = Mesh.Elements(i,2);
            Mesh.Elements(i,2) = tmp;
        end
    end

    % Plot shap function

    OFFSET = 0.05;

    % Compute axes limits

    XMin = min(Mesh.Coordinates(:,1));
    XMax = max(Mesh.Coordinates(:,1));
    YMin = min(Mesh.Coordinates(:,2));
    YMax = max(Mesh.Coordinates(:,2));
    XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
    YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];

    % Compute color axes limits

    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)          % or error will occur in set function
        CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
        CLim = [1-OFFSET 1+OFFSET]*CMin;
    end

    % Plot figure
    
    set(gcf,'Color',[1 1 1]*p);
    set(gca,'CameraPosition', [1 -.5 1.3]*50,...
        'CameraTarget',[0 0 0], ...
        'CameraUpVector',[0 0 1], ...
        'CameraViewAngle',1.4, ...
        'DataAspectRatio', [1 1 0.9],...
        'Position',[0 0 1 1], ...
        'Visible','on', ...
            'XLim',[-0.1 1.1], ...
            'YLim',[-0.1 1.1], ...
            'ZLim',[CMin,CMax]);
         xlabel('\bfX')
        ylabel('\bfY')
        zlabel('\bfZ')
    patch('faces', Mesh.Elements, ...
        'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
        'CData', U, ...
        'EdgeColor',[0,0,0], ...
        'FaceColor','flat', ...
        'FaceLighting','phong', ...
        'AmbientStrength',0.3, ...
        'DiffuseStrength',0.6, ...
        'Clipping','off',...
        'BackFaceLighting','lit', ...
        'SpecularStrength',1.1, ...
        'SpecularColorReflectance',1, ...
        'SpecularExponent',7);
    alpha(.7)
   % l1 = light('Position',[-1 -.5 1.3], ...
%             'Style','infinite', ...
%             'Color',[0 0.8 0.8]);
%       l2 = light('Position',[1 1 2], ...
%             'Style','local', ...
%             'Color',[0.8 .8 0]);
%         l3 = light('Position',[1 -1 .5], ...
%             'Style','local', ...
%             'Color',[0.8 .8 0]);
    
    if(nargout > 0)
        varargout{1} = fig;
    end
    
return


    