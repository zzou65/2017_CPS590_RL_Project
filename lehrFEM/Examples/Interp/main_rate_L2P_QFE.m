% Run script for interest

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constant
    
    NREFS = 3;
    lambda = linspace(1.9,2.1,5);
    lambda = [0.1:0.2:1.7 lambda 2.2:0.2:3.0];
    NLAMBDA = length(lambda);
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,NLAMBDA);
    Err_H1 = zeros(NREFS,NLAMBDA);
    p = zeros(NLAMBDA,2);
    
    % Initialize mesh

    X0 = [0 0];                                     % Lower left corner point of rectangle
    A = 1;                                          % Length of rectangle
    BBOX = [X0; X0+[A A]];                          % Bounding box
    H0 = 0.2;                                       % Initial mesh width
    DHANDLE = @dist_rect;                           % Signed distance function
    HHANDLE = @h_uniform;                           % Element size function
    FIXEDPOS = [X0; X0+[A 0]; X0+[A A]; X0+[0 A]];  % Fixed boundary vertices of the mesh
    DISP = 0;                                       % Display flag
    Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,X0,A,A);
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
    
    % Run the script
    
    for iter = 1:NREFS
        
        % Do REG refine
        
        Mesh = refine_REG(Mesh);
        
        % Get mesh width
        
        M_W(iter) = get_MeshWidth(Mesh);
        
        for j = 1:NLAMBDA

            % Obtain singular Funtion

            F_HANDLE = @(x,varargin)(x(:,1).^2+x(:,2).^2).^(lambda(j)/2);
            GD_HANDLE = @(x,varargin)[lambda(j)*x(:,1).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1) ...
                                      lambda(j)*x(:,2).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1)];
            
            
            % Compute the energy error

            U = L2_interp_QFE(Mesh,P7O6(),F_HANDLE);
            Err_L2(iter,j) = L2Err_QFE(Mesh,U,P7O6(),F_HANDLE);
            Err_H1(iter,j) = H1SErr_QFE(Mesh,U,P7O6(),GD_HANDLE);
                        
        end
        
    end
    
    % Compute slope
    
    for i = 1:NLAMBDA
        
        slope = polyfit(log(M_W),log(Err_L2(:,i)),1);
        p(i,1) = slope(1);
        slope = polyfit(log(M_W),log(Err_H1(:,i)),1);
        p(i,2) = slope(1);
        
    end
    p(find(p<0))=NaN;
    
    % Generate .eps files
    
    figure;
    plot(lambda,p(:,1),'-x',...
         lambda,p(:,2),'->'); 
    grid on
    title('\bfRelationship between convergence rate and lambda(L2P QFE)');
    xlabel('\bf\lambda');
    ylabel('\bfConvergence rate');
    axis([0,3.5,0,3.5])
    legend('\bfL2','\bfH1',...
            'Location','NorthEast')
%     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
%     add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_L2P_QFE.eps');
    !gv rate_L2P_QFE.eps &
    
    % Clear memory
    
%     clear all;