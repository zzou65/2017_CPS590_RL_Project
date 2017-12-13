% Run script for interest

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constant
    
    NREFS = 7;
    lambda = 0.2:0.2:2.8;
    NLAMBDA = length(lambda);
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,NLAMBDA);
    Err_H1 = zeros(NREFS,NLAMBDA);
    p = zeros(NLAMBDA,2);
    
    save data_LFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA
    
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
%          figure;
%          plot_Mesh(Mesh,'as');
%          file = ['MeshLFE_' int2str(iter) '.eps'];
%          print('-depsc',file);
%    
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
            QuadRule = Duffy(TProd(gauleg(0,1,10)));
            U = LFE_DOF_interp(Mesh,F_HANDLE);
            Err_L2(iter,j) = L2Err_LFE(Mesh,U,QuadRule,F_HANDLE);
            Err_H1(iter,j) = H1SErr_LFE(Mesh,U,QuadRule,GD_HANDLE);
            clear U;
                        
        end
        save data_LFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA -APPEND
    end
    % Generate .eps files
    
    save data_LFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA -APPEND
    % L2 error
    figure;
    plot(M_W,Err_L2(:,2),'-x',...
         M_W,Err_L2(:,3),'-^',...
         M_W,Err_L2(:,4),'-h',...
         M_W,Err_L2(:,5),'-o',...
         M_W,Err_L2(:,10),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bfh','FontSize',14);
    ylabel('\bfinterp. error (L^2 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(2))],['\bf\lambda=' num2str(lambda(3))],...
           ['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(5))],...
           ['\bf\lambda=' num2str(lambda(10))],...
            'Location','NorthEast')
    %p = polyfit(log(M_W),log(Err_L2(:,1)),1);
    %add_Slope(gca,'SouthEast',2);
    print('-depsc', 'rate_L2_LFE.eps');
    % !gv rate_L2_LFE-DOF.eps &
    
    
    % H1 error
    figure;
    plot(M_W,Err_H1(:,2),'-x',...
         M_W,Err_H1(:,3),'-^',...
         M_W,Err_H1(:,4),'-h',...
         M_W,Err_H1(:,5),'-o',...
         M_W,Err_H1(:,10),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bfh','FontSize',14);
    ylabel('\bfinterp. error (H^1 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(2))],['\bf\lambda=' num2str(lambda(3))],...
           ['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(5))],...
           ['\bf\lambda=' num2str(lambda(10))],...
            'Location','NorthEast')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_H1_LFE.eps');
    %!gv rate_H1_LFE-DOF.eps &
    
   
    % Estimate for Convergence rate
    
    est_L2=zeros((NREFS-1),5);
    est_H1=zeros((NREFS-1),5);
    for i=2:NREFS
       est_L2(i-1,:)=log(Err_L2(i,[2,3,4,5,10])./Err_L2(i-1,[2,3,4,5,10]))./log(M_W(i)./M_W(i-1));
       est_H1(i-1,:)=log(Err_H1(i,[2,3,4,5,10])./Err_H1(i-1,[2,3,4,5,10]))./log(M_W(i)./M_W(i-1));
    end
    figure;
    plot(2:NREFS,est_L2(:,1),'x',...
         2:NREFS,est_L2(:,2),'^',...
         2:NREFS,est_L2(:,3),'h',...
         2:NREFS,est_L2(:,4),'o',...
         2:NREFS,est_L2(:,5),'s'); 
    grid on
    %set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bf level','FontSize',14);
    ylabel('\bf est.rate (L^2 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(2))],['\bf\lambda=' num2str(lambda(3))],...
           ['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(5))],...
           ['\bf\lambda=' num2str(lambda(10))],...
            'Location','NorthWest')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_L2_est_LFE.eps');
    %!gv rate_H1_LFE-DOF.eps &
    figure;
    plot(2:NREFS,est_H1(:,1),'x',...
         2:NREFS,est_H1(:,2),'^',...
         2:NREFS,est_H1(:,3),'h',...
         2:NREFS,est_H1(:,4),'o',...
         2:NREFS,est_H1(:,5),'s'); 
    grid on
    %set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bf level','FontSize',14);
    ylabel('\bf est. rate (H^1 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(2))],['\bf\lambda=' num2str(lambda(3))],...
           ['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(5))],...
           ['\bf\lambda=' num2str(lambda(10))],...
            'Location','NorthWest')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_H1_est_LFE.eps');
    %!gv rate_H1_LFE-DOF.eps &
    
   
    
    % CONVERGENCE RATE
    
    % Compute slope
    
    for i = 1:NLAMBDA
        
        slope = polyfit(log(M_W((NREFS-2):NREFS)),log(Err_L2(((NREFS-2):NREFS),i)),1);
        p(i,1) = slope(1);
        slope = polyfit(log(M_W((NREFS-2):NREFS)),log(Err_H1(((NREFS-2):NREFS),i)),1);
        p(i,2) = slope(1);
        
    end
    
    figure;
    plot(lambda,p(:,1),'-x',...
         lambda,p(:,2),'-^'); 
    grid on
    xlabel('\bf\lambda','FontSize',14);
    ylabel('\bfconvergence rate','FontSize',14);
    axis([0,3,0,2.5])
    legend('\bfL^2','\bfH^1',...
            'Location','NorthEast')
%     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
%     add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_LFE.eps');
    !gv rate_H1P_LFE.eps &
    
    % Clear memory
    clear all;
