% Run script for all the six QFE shap functions with transparent effect and 
% without outlines 

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Intialize constant
    
    NMesh = 100;
    TOL = 1e-10;
    
    % Generate data
    
    x = linspace(0,1,NMesh);
    y = linspace(0,1,NMesh);
    [X,Y] = meshgrid(x,y);
    loc = find(X+Y<=1);
    X = X(loc);
    Y = Y(loc);
    shap = shap_LFE([X Y]);
    
    for i=1:3

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap(loc,i));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap(loc,i)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap(:,i));
        
        % Output .eps files

        FileName = ['Shap_LFE3D_' int2str(i) '.eps'];
        print('-depsc', FileName);
        system(['gv ' FileName ' &']);
        
    end
    
    % Clear Memory
    
    clear all;