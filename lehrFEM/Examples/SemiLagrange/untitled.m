for i=1:10
    nsteps=2^i;
    y=zeros(nsteps+1,1);
    t=zeros(nsteps+1,1);
    y(1)=1;
    t(1)=0;
    h=5/nsteps;
    for k=1:(nsteps)
        t(k+1)=t(k)+h;
        y(k+1)=y(k)+h*(-1/(1+t(k+1))^2);
    end
    plot(t,abs(y-1./(1+t)));
    hold on;
    err(i)=abs(y(end)-1./(1+t(end)));
    s(i)=h;
end
set(gca,'YScale','log');
hold off;
figure;
plot(s,err);
grid('on');
set(gca,'YScale','log','XScale','log');
p = polyfit(log(s),log(err),1);
add_Slope(gca,'East',p(1));
    