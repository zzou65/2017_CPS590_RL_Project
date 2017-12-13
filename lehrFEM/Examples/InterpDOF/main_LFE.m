% Run script for error of H1 projection using linear finite element

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumannn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    NREFS = 7;
    lambda = [0.4 0.5 0.8 1.0 2.0];
    NLAMBDA = length(lambda);
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,NLAMBDA);
    Err_H1 = zeros(NREFS,NLAMBDA);
    
    % Initialize mesh

    X0 = [-1 -1];                                     % Lower left corner point of rectangle
    A = 2;                                          % Length of rectangle
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
        
        Mesh = refine_REG_jiggle(Mesh);
        %Mesh = refine_REG(Mesh);
        % Get mesh width
        
        M_W(iter) = get_MeshWidth(Mesh);
        
        for j = 1:NLAMBDA

            % Obtain singular Funtion

            F_HANDLE = @(x,varargin)(x(:,1).^2+x(:,2).^2).^(lambda(j)/2);
            GD_HANDLE = @(x,varargin)[lambda(j)*x(:,1).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1) ...
                                      lambda(j)*x(:,2).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1)];
            
            
            % Compute the energy error

            U = LFE_DOF_interp(Mesh,F_HANDLE);
            Err_L2(iter,j) = L2Err_LFE(Mesh,U,P7O6(),F_HANDLE);
            Err_H1(iter,j) = H1SErr_LFE(Mesh,U,P7O6(),GD_HANDLE);

        end
        
    end
    
    % Generate .eps files
    
    figure;
    plot(M_W,Err_L2(:,1),'-x',...
         M_W,Err_L2(:,2),'-^',...
         M_W,Err_L2(:,3),'-h',...
         M_W,Err_L2(:,4),'-o',...
         M_W,Err_L2(:,5),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for L2 error of LFE-DOFs interpolation');
    xlabel('\bfh','fontsize=14pt');
    ylabel('\bfEnergy error','fontsize=4pt');
    legend(['\bf\lambda=' num2str(lambda(1))],['\bf\lambda=' num2str(lambda(2))],...
           ['\bf\lambda=' num2str(lambda(3))],['\bf\lambda=' num2str(lambda(4))],...
           ['\bf\lambda=' num2str(lambda(5))],...
            'Location','NorthEast')
%     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
    %add_Slope(gca,'SouthEast',2);
    print('-depsc', 'rate_L2_LFE.eps');
   % !gv rate_L2_LFE-DOF.eps &
    
    figure;
    plot(M_W,Err_H1(:,1),'-x',...
         M_W,Err_H1(:,2),'-^',...
         M_W,Err_H1(:,3),'-h',...
         M_W,Err_H1(:,4),'-o',...
         M_W,Err_H1(:,5),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for H1 error of LFE-DOFs interpolation');
    xlabel('\bfh','fontsize=14pt');
    ylabel('\bfEnergy error','fontsize=14pt');
    legend(['\bf\lambda=' num2str(lambda(1))],['\bf\lambda=' num2str(lambda(2))],...
           ['\bf\lambda=' num2str(lambda(3))],['\bf\lambda=' num2str(lambda(4))],...
           ['\bf\lambda=' num2str(lambda(5))],...
            'Location','NorthEast')
%     p = polyfit(log(M_W),log(Err_H1(:,1)),1);
   % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_H1_LFE.eps');
    %!gv rate_H1_LFE-DOF.eps &
    
    
    % Clear memory
    
    clear all;
    
    
    