% Run script for all the six QFE shap functions with transparent effect and 
% without outlines 

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Intialize constant
    
    NMesh = 25;
    TOL = 1e-10;
    
    % Generate data
    
    x = linspace(0,1,NMesh);
    y = linspace(0,1,NMesh);
    [X,Y] = meshgrid(x,y);
    loc = find(X+Y<=2);
    X = X(loc);
    Y = Y(loc);
    shap = shap_BFE([X Y]);
    
    for i=1:4
        % billinar mapping from unit square to unit triangle
        X_t=X-0.5.*(X.*Y);
        Y_t=Y-0.5.*(X.*Y);
        % Plot figures
%         figure
%         plot3(X_t,Y_t,shap(:,i));
%         grid on;
        figure
        plot([0 1 0 0],[0 0 1 0],'-rs','linewidth',1.5,'MarkerEdgeColor',...
            'k','MarkerFaceColor','r','MarkerSize',10)
        hold on
        X_t=X-0.5.*(X.*Y);
        Y_t=Y-0.5.*(X.*Y);
        loc = find(sum([X_t Y_t].^2,2) < TOL | sum([X_t-1 Y_t].^2,2) < TOL | sum([X_t Y_t-1].^2,2) < TOL | sum([X_t-0.5 Y_t-0.5].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap(loc,i));
        loc = loc(pos);
        plot3([X_t(loc) X_t(loc)] ,[Y_t(loc) Y_t(loc)],[zeros(size(loc)) shap(loc,i)],'--k','linewidth',1.5)
        plot_Shap([X_t Y_t],shap(:,i));
        
        % Output .eps files

%       FileName = ['Shap_LFE3D_' int2str(i) '.eps'];
%        print('-depsc', FileName);
%        system(['gv ' FileName ' &']);
        
    end
    
    % Clear Memory
    
    clear all;