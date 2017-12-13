  load data2
  E2=Err2;
  i2=iter2;
  n2=nDofs2;
  load data;
  Err2=E2;
  iter2=i2;
  nDofs2=n2;
  % Generate figures
    
  fig = figure('Name','Recovery based error estimator');
  plot(nDofs,Err,'b-*',nDofs1,Err1,'r-o',nDofs2,Err2,'g-^');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs }');
  ylabel('{\bf Errors}');
  legend('xxxxxxxxxxxxxxx','yyyyyyyyyyyyyyy','zzzzzzzzzzzzzzz','Location','SouthWest');
  p = polyfit(log(nDofs((iter-10):iter)),log(Err((iter-10):iter)),1);
  add_Slope(gca,'South',p(1));
  p = polyfit(log(nDofs1((iter1-10):iter1)),log(Err1((iter1-10):iter1)),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(nDofs2((iter2-4):iter2)),log(Err2((iter2-4):iter2)),1);
  add_Slope(gca,'East',p(1));
   print('-depsc','adap.eps');
  % Clear memory