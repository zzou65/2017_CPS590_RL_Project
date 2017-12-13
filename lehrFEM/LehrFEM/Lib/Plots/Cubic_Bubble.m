function varargout = Cubic_Bubble(Vertices)
% CUBIC_BUBBLE Plot cubic bubble.
%
%   CUBIC_BUBBLE(VERTICES) generates a plot of cubbic bubble on the region
%   represented by Vertices.
%
%   H = CUBIC_BUBBLE(VERTICES) also returns the handle to the figure.
%
%   Example:
%
%   Cubic_Bubble([0 0;1 0;0 1]);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    NRefs = 4;
    BK = [Vertices(2,:)-Vertices(1,:);Vertices(3,:)-Vertices(1,:)];
    det_BK = abs(det(BK));
  
    Mesh.Coordinates = Vertices;
    Mesh.Elements = [1 2 3];

    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
  
    for i = 1:NRefs

        Mesh = refine_REG(Mesh);

    end

    x1 = Vertices(1,1);  y1 = Vertices(1,2);
    x2 = Vertices(2,1);  y2 = Vertices(2,2);
    x3 = Vertices(3,1);  y3 = Vertices(3,2);



    F_Handle = @(x,varargin)27/det_BK^3* ...
        ((x(:,1)-x2)*(y2-y3)+(x(:,2)-y2)*(x3-x2)).* ...
        ((x(:,1)-x3)*(y3-y1)+(x(:,2)-y3)*(x1-x3)).* ...
        ((x(:,1)-x1)*(y1-y2)+(x(:,2)-y1)*(x2-x1));


    M = assemMat_LFE(Mesh,@MASS_LFE);
    L = assemLoad_LFE(Mesh,P7O6,F_Handle);
    U = M\L;

    fig = plot([x1 x2 x3 x1],[y1 y2 y3 y1],'.-r','LineWidth',1.5)
    set(gca,'CameraPosition', [x3+(x1-x2) y3+(y1-y2) 1.6],...
        'CameraTarget',[(x1+x2+x3)/3 (y1+y2+y3)/3 0], ...
        'Visible','off')
    h = patch('faces', Mesh.Elements, ...
        'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
        'CData', U, ...
        'facecolor', [.1 0 0], ...
        'edgecolor', [1 0 0],...
        'LineWidth',1.5);
    alpha(.3)

    if(nargout > 0)
        varargout{1} = fig;
    end

    
return
    
    