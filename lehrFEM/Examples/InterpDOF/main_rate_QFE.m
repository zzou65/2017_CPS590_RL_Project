% Run script for interest

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constant
    
    NREFS =7;
    lambda = 0.2:0.2:2.8;
    %lambda = linspace(1.9,2.1,5);
    lambda = [0.2:0.2:1.8 1.9 2.2:0.2:2.8];
    NLAMBDA = length(lambda);
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,NLAMBDA);
    Err_H1 = zeros(NREFS,NLAMBDA);
    p = zeros(NLAMBDA,2);
    save data_QFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA
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
%         figure;
%         plot_Mesh(Mesh,'as');
%         file = ['MeshQFE_' int2str(iter) '.eps'];
%         print('-depsc',file);
        
        % Get mesh width
        
        M_W(iter) = get_MeshWidth(Mesh);
        
        for j = 1:NLAMBDA

            % Obtain singular Funtion

            F_HANDLE = @(x,varargin)(x(:,1).^2+x(:,2).^2).^(lambda(j)/2);
            GD_HANDLE = @(x,varargin)[lambda(j)*x(:,1).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1) ...
                                      lambda(j)*x(:,2).*(x(:,1).^2+x(:,2).^2).^(lambda(j)/2-1)];
            
            
            % Compute the energy error
            QuadRule = Duffy(TProd(gauleg(0,1,10)));
            U = QFE_DOF_interp(Mesh,F_HANDLE);
            Err_L2(iter,j) = L2Err_QFE(Mesh,U,QuadRule,F_HANDLE);
            Err_H1(iter,j) = H1SErr_QFE(Mesh,U,QuadRule,GD_HANDLE);
            clear U;            
        end
       
        save data_QFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA -APPEND
    end
    save data_QFE M_W Err_L2 Err_H1 lambda NREFS NLAMBDA -APPEND
    % Generate .eps files
    
    % L2 Norm error
    figure;
    plot(M_W,Err_L2(:,4),'-x',...
         M_W,Err_L2(:,6),'-^',...
         M_W,Err_L2(:,8),'-h',...
         M_W,Err_L2(:,10),'-o',...
         M_W,Err_L2(:,12),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bfh','FontSize',14);
    ylabel('\bfinterp. error (L^2 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(6))],...
           ['\bf\lambda=' num2str(lambda(8))],['\bf\lambda=' num2str(lambda(10))],...
           ['\bf\lambda=' num2str(lambda(12))],...
            'Location','NorthEast')
    %     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
    %add_Slope(gca,'SouthEast',2);
    print('-depsc', 'rate_L2_QFE.eps'); 
    %    !gv rate_L2_QFE-DOF.eps &
    
    % H1 norm error
    figure;
    plot(M_W,Err_H1(:,4),'-x',...
         M_W,Err_H1(:,6),'-^',...
         M_W,Err_H1(:,8),'-h',...
         M_W,Err_H1(:,10),'-o',...
         M_W,Err_H1(:,12),'-s'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    xlabel('\bfh','FontSize',14);
    ylabel('\bfinterp. error (H^1 norm)','FontSize',14);
    legend(['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(6))],...
           ['\bf\lambda=' num2str(lambda(8))],['\bf\lambda=' num2str(lambda(10))],...
           ['\bf\lambda=' num2str(lambda(12))],...
            'Location','NorthEast')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    %add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_H1_QFE.eps');
    %!gv rate_H1_QFE-DOF.eps &
    
    
    
    % Estimate for Convergence rate
    
    est_L2=zeros((NREFS-1),5);
    est_H1=zeros((NREFS-1),5);
    for i=2:NREFS
       est_L2(i-1,:)=log(Err_L2(i,[4,6,8,10,12])./Err_L2(i-1,[4,6,8,10,12]))./log(M_W(i)./M_W(i-1));
       est_H1(i-1,:)=log(Err_H1(i,[4,6,8,10,12])./Err_H1(i-1,[4,6,8,10,12]))./log(M_W(i)./M_W(i-1));
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
    legend(['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(6))],...
           ['\bf\lambda=' num2str(lambda(8))],['\bf\lambda=' num2str(lambda(10))],...
           ['\bf\lambda=' num2str(lambda(12))],...
            'Location','NorthWest')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_L2_est_QFE.eps');
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
    legend(['\bf\lambda=' num2str(lambda(4))],['\bf\lambda=' num2str(lambda(6))],...
           ['\bf\lambda=' num2str(lambda(8))],['\bf\lambda=' num2str(lambda(10))],...
           ['\bf\lambda=' num2str(lambda(12))],...
            'Location','NorthWest')
    %p = polyfit(log(M_W),log(Err_H1(:,1)),1);
    % add_Slope(gca,'SouthEast',1);
    print('-depsc', 'rate_H1_est_QFE.eps');
    %!gv rate_H1_LFE-DOF.eps &
    
  
    % Convergence rate
    
    % Compute slope
    
    for i = 1:NLAMBDA
        
        slope = polyfit(log(M_W((NREFS-2):NREFS)),log(Err_L2(((NREFS-2):NREFS),i)),1);
        p(i,1) = slope(1);
        slope = polyfit(log(M_W((NREFS-2):NREFS)),log(Err_H1(((NREFS-2):NREFS),i)),1);
        p(i,2) = slope(1);
        
    end
    p(find(p<0))=NaN;
        
    
    figure;
    plot(lambda,p(:,1),'-*',...
         lambda,p(:,2),'->'); 
    grid on
    xlabel('\bf\lambda','FontSize',14);
    ylabel('\bfconvergence rate','FontSize',14);
    axis([0,3.5,0,3])
    legend('\bfL^2','\bfH^1',...
            'Location','NorthEast')
%     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
%     add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_QFE.eps');
    !gv rate_H1P_QFE.eps &
    
    % Clear memory
    
    clear all;