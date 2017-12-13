%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u1 sq   %%%%%%%%%%%%%%%%%%%

M_W=[.3536 .1768 .0884 .0442 .0221];

load MIXDG11.mat L2Err
Err_L2_1=L2Err;
load WReg11.mat L2Err
Err_L2_3=L2Err;
load SReg11.mat L2Err
Err_L2_4=L2Err;
figure
    
    plot(M_W,Err_L2_1,'-^',M_W,Err_L2_3,'-o',M_W,Err_L2_4,'-+');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'South',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_3)),1);
    add_Slope(gca,'North',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_4)),1);
    add_Slope(gca,'West',p(1));
    ylabel('\bf\fontsize{14} Relative discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for relative discretization error on the square domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F','\bf\fontsize{10}SRegLFE')
   %  Generate .eps files
    
    print('-depsc', 'rate_u1sq.eps');
    !gv rate_u1sq.eps &
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u2 ls   %%%%%%%%%%%%%%%%%%%%    
load MIXDG23.mat L2Err
Err_L2_1=L2Err;
load WReg23.mat L2Err
Err_L2_3=L2Err;
load SReg23.mat L2Err
Err_L2_4=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-^',M_W,Err_L2_3,'-o',M_W,Err_L2_4,'-+');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'South',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_3)),1);
    add_Slope(gca,'West',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_4)),1);
    add_Slope(gca,'North',p(1));
    ylabel('\bf\fontsize{14} Relative discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for relative discretization error on the L-shaped domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F','\bf\fontsize{10}SRegLFE')
   %  Generate .eps files
    
    print('-depsc', 'rate_u2ls.eps');
    !gv rate_u2ls.eps &

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    u3 ls   %%%%%%%%%%%%%%%%%%%%
load MIXDG24.mat L2Err
Err_L2_1=L2Err;
load WReg24.mat L2Err
Err_L2_3=L2Err;
load SReg24.mat L2Err
Err_L2_4=L2Err;

figure
    
    plot(M_W,Err_L2_1,'-^',M_W,Err_L2_3,'-o',M_W,Err_L2_4,'-+');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2_1)),1);
    add_Slope(gca,'South',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_3)),1);
    add_Slope(gca,'West',p(1));
    p = polyfit(log(M_W),log(abs(Err_L2_4)),1);
    add_Slope(gca,'NorthEast',p(1));
    ylabel('\bf\fontsize{14} Relative discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for relative discretization error on the L-shaped domain')
    legend('\bf\fontsize{10}MIXDG','\bf\fontsize{10}WRegW1F','\bf\fontsize{10}SRegLFE')
    
   %  Generate .eps files
    
    print('-depsc', 'rate_u3ls.eps');
    !gv rate_u3ls.eps &
    