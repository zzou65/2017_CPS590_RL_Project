function example5()

load ./results/BITexample5_0.8.mat
err=stopTimeL2ErrH(:,1);
load ./results/BITexample5_0.6.mat
err=[err,stopTimeL2ErrH(:,1)];
load ./results/BITexample5_0.4.mat
err=[err,stopTimeL2ErrH(:,1)];
load ./results/BITexample5_0.2.mat
err=[err,stopTimeL2ErrH(:,1)];

figure; 
hold on;
plot(mw,err(:,1),'x-',...
    mw,err(:,2),'x-',...
    mw,err(:,3),'x-',...
    mw,err(:,4),'x-',...
    mw,err(1,1)/mw(1)*mw,'k--','Linewidth',2,'MarkerSize',8);
grid('on');
set(gca,'YScale','log','XScale','log','FontSize',12,'FontWeight','bold');
xlabel('{\bf h}','FontSize',12,'FontWeight','bold');
ylabel('{\bf }','FontSize',12,'FontWeight','bold');
legend('CFL=0.8','CFL=0.6','CFL=0.4','CFL=0.2','h','Location','Northwest');
hold off;

saveas(gcf,'./results/BITexample5.fig')
print('-depsc' ,'./results/BITexample5.eps')

