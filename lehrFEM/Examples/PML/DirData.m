function gd=DirData(x,Flag)
   if (Flag==-1) gd=besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
       %gd=exp(i*sqrt(x(:,1).^2+x(:,2).^2));
   else %gd=besselh(0,sqrt(x(:,1).^2+x(:,2).^2));
         gd=zeros(size(x,1),1);
   end
  return