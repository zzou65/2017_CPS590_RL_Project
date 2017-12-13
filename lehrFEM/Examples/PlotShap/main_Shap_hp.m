% Run script for all hp shap functions with transparent effect 
% without outlines 

%   Copyright 2005-2008 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Intialize constant
    
    NMesh = 100;
    TOL = 1e-10;
    pmax=3;
    
    % Generate data
    
    x = linspace(0,1,NMesh);
    y = linspace(0,1,NMesh);
    [X,Y] = meshgrid(x,y);
    loc = find(X+Y<=1);
    X = X(loc);
    Y = Y(loc);
    shap = shap_hp([X Y],pmax);
    
    % vertex
    for i=1:3

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap.vshap.shap{i}(loc));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap.vshap.shap{i}(loc)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap.vshap.shap{i});
            
    end
    
    % edge 1 
    Edofs1=size(shap.eshap{1}.shap,2)
    
    for i=1:Edofs1

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap.eshap{1}.shap{i}(loc));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap.eshap{1}.shap{i}(loc)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap.eshap{1}.shap{i});
            
    end
    
    % edge 2 
    Edofs2=size(shap.eshap{2}.shap,2)
    
    for i=1:Edofs2

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap.eshap{2}.shap{i}(loc));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap.eshap{2}.shap{i}(loc)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap.eshap{2}.shap{i});
            
    end
    
    % edge 3 
    Edofs3=size(shap.eshap{3}.shap,2)
    
    for i=1:Edofs3

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap.eshap{3}.shap{i}(loc));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap.eshap{3}.shap{i}(loc)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap.eshap{3}.shap{i});
            
    end
    
     % element 
    Cdofs=size(shap.cshap.shap,2)
    
    for i=1:Cdofs

        % Plot figures
        
        figure
        plot([0 1 0 0],[0 0 1 0],'.-k','linewidth',1.5,'markersize',10)
        hold on
        loc = find(sum([X Y].^2,2) < TOL | sum([X-1 Y].^2,2) < TOL | sum([X Y-1].^2,2) < TOL);
        [pos dummy1 dummy2] = find(shap.cshap.shap{i}(loc));
        loc = loc(pos);
        plot3([X(loc) X(loc)] ,[Y(loc) Y(loc)],[zeros(size(loc)) shap.cshap.shap{i}(loc)],'--k','linewidth',1.5)
        plot_Shap([X Y],shap.cshap.shap{i});
            
    end
    
    % Clear Memory
    
    clear all;