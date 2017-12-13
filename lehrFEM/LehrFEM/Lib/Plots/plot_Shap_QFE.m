function varargout = plot_Shap_QFE(NMesh,QNumber,Vertex)
% PLOT_SHAP_QFE plot QFE shape functions in 3D with lighting effect
%
%   PLOT_SHAP_QFE(NMESH,LNUMBER,VERTEX) generates a 3D plot of the QFE
%   shape function on the domain represented by Vertex
%
%   The number NMesh determines how fine the mesh is.
%   
%   The number LNumber, which is integer from 1 to 6, determines which QFE
%   shape function to take.
%   
%   The matrix Vertex determine the set of point which represents the
%   destination of the mapping.
%
%   H = PLOT_SHAP_QFE(NMESH,LNUMBER,VERTEX) also returns the handle to the 
%   figure.
%
%   Example:
%
%   plot_Shap_QFE(100,1,[0 0;1 0;0 1]);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Intialize constant
    
    TOL = 1e-10;
    
    % Generate data on reference triangle domain
    
    x = linspace(0,1,NMesh);
    y = linspace(0,1,NMesh);
    [X,Y] = meshgrid(x,y);
    loc = find(X+Y<=1);
    X = X(loc);
    Y = Y(loc);
    shap = shap_QFE([X Y]);
    NV = affine_map([X Y],Vertex);
    X = NV(:,1);
    Y = NV(:,2);
    
    % Plot outline
    
    px1 = linspace(0,1,NMesh)';
    py1 = zeros(size(px1));
    pz1 = shap_QFE([px1 py1]);
    py2 = linspace(0,1,NMesh)';
    px2 = zeros(size(py2));
    pz2 = shap_QFE([px2 py2]);
    px3 = linspace(0,1,NMesh)';
    py3 = linspace(1,0,NMesh)';
    pz3 = shap_QFE([px3 py3]);
    V1 = affine_map([px1 py1],Vertex);
    V2 = affine_map([px2 py2],Vertex);
    V3 = affine_map([px3 py3],Vertex);
    plot3(V1(:,1),V1(:,2),pz1(:,QNumber),'-k','linewidth',1.5);
    hold on
    plot3(V2(:,1),V2(:,2),pz2(:,QNumber),'-k','linewidth',1.5);
    hold on
    plot3(V3(:,1),V3(:,2),pz3(:,QNumber),'-k','linewidth',1.5);
    hold on
    
    % Plot figure

    plot([Vertex(:,1);Vertex(1,1)],[Vertex(:,2);Vertex(1,2)],'.-k','linewidth',1.5,'markersize',10)
    hold on
    loc = find(sum([X-Vertex(1,1) Y-Vertex(1,2)].^2,2) < TOL |...
               sum([X-Vertex(2,1) Y-Vertex(2,2)].^2,2) < TOL |...
               sum([X-Vertex(3,1) Y-Vertex(3,2)].^2,2) < TOL);
    [pos dummy1 dummy2] = find(shap(loc,QNumber));
    loc = loc(pos);
    plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap(loc,QNumber)],'--k','linewidth',2)
    plot_Shap([X Y],shap(:,QNumber));

    if(nargout > 0)
        varargout{1} = fig;
    end
    
return