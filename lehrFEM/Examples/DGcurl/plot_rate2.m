%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u1 sq   %%%%%%%%%%%%%%%%%%%

M_W=[.3536 .1768 .0884 .0442 .0221];

load DGcurl11.mat L2Err
Err_L2_1=L2Err;
load WReg11.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} L2 discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the square domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u1sq.eps');
    !gv rate_u1sq.eps &
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u2 sq   %%%%%%%%%%%%%%%%%%%  
load DGcurl13.mat L2Err
Err_L2_1=L2Err;
load WReg13.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'NorthEast',p(1));
    ylabel('\bf\fontsize{14} L2 discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the square domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u2sq.eps');
    !gv rate_u2sq.eps &
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u3 sq   %%%%%%%%%%%%%%%%%%%
load DGcurl14.mat L2Err
Err_L2_1=L2Err;
load WReg14.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} L2 discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the square domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u3sq.eps');
    !gv rate_u3sq.eps &
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u1 ls   %%%%%%%%%%%%%%%%%%%%
M_W=[.5 .25 .125 .0625 .03125];

load DGcurl21.mat L2Err
Err_L2_1=L2Err;
load WReg21.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} L2 discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the L-shaped domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u1ls.eps');
    !gv rate_u1ls.eps &
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u2 ls   %%%%%%%%%%%%%%%%%%%%    
load DGcurl23.mat L2Err
Err_L2_1=L2Err;
load WReg23.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} L2  discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the L-shaped domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u2ls.eps');
    !gv rate_u2ls.eps &

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u3 ls   %%%%%%%%%%%%%%%%%%%%
load DGcurl24.mat L2Err
Err_L2_1=L2Err;
load WReg24.mat L2Err
Err_L2_2=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-o',M_W,Err_L2_2,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'SouthEast',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_2)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} L2 discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for L2 discretization error on the L-shaped domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F')
   %  Generate .eps files
    
    print('-depsc', 'rate_u3ls.eps');
    !gv rate_u3ls.eps &
    