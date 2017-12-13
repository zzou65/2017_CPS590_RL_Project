% plot_cornerwave.m                 Andrea Moiola           5.5.2009
% script to plot the exact solution
% J_xi(omega r) * cos(xi*theta)
% in the right unit square, surf o contour

clear all;
np=100;
R=1;
x=linspace(0,R,np);
y=linspace(-R/2,R/2,np);
[xx,yy]=meshgrid(x,y);

xi=3/2
omega=10;
theta_0=0;
lambda=omega;
vangle = @(x) angle(x(:,1)+sqrt(-1)*x(:,2));
vnorm=@(x) sqrt(x(:,1).^2+x(:,2).^2) ;          %euclid norm for array of 2D vector
u_ex = @(x,varargin) (besselj(xi, lambda*vnorm(x) )) .*cos(xi*(vangle(x)-theta_0));


zz=zeros(np,np);
for j=1:np
    for k=1:np
        zz(j,k)=u_ex([xx(j,k), yy(j,k)]);
    end
end

figure
%surf(xx,yy,zz)
%shading interp
% %shading flat

contour(xx,yy,zz, 'linewidth', 2)

view([0,0,1])
colorbar
axis equal
set(gca, 'FontSize', 18);
