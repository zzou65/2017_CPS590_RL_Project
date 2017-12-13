% Run script for interest

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

load ~/local/interpol/data_QFE.mat
%NREFS=6;
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
     plot(1:(NREFS-1),est_L2(:,1),'x',...
         1:(NREFS-1),est_L2(:,2),'^',...
         1:(NREFS-1),est_L2(:,3),'h',...
         1:(NREFS-1),est_L2(:,4),'o',...
         1:(NREFS-1),est_L2(:,5),'s'); 
    grid on
    set(gca,'XTick',[0,1:6,7],'XTickLabel',{' ','1','2','3','4','5','6',''},'XLim',[0.5,6.5])
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
    plot(1:(NREFS-1),est_H1(:,1),'x',...
         1:(NREFS-1),est_H1(:,2),'^',...
         1:(NREFS-1),est_H1(:,3),'h',...
         1:(NREFS-1),est_H1(:,4),'o',...
         1:(NREFS-1),est_H1(:,5),'s'); 
    grid on
    set(gca,'XTick',[0,1:6,7],'XTickLabel',{' ','1','2','3','4','5','6',''},'XLim',[0.5,6.5])
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
    axis([0,3,0,3.5])
    legend('\bfL^2','\bfH^1',...
            'Location','SouthEast')
%     p = polyfit(log(M_W),log(Err_L2(:,1)),1);
%     add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_QFE.eps');
    !gv rate_H1P_QFE.eps &
    
    % Clear memory
    
    clear all;