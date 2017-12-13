function e_f(n,m)
% E_F(N,M) plots eigen functions
%
%   E_F(N,M) plots all the four eigen functions with Nth in x direction and 
%   Mth in y deriction on the square domain
%   
%   Example:
%
%   e_f(1,2);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    figure
    
    x = linspace(-1,1,100);
    y = linspace(-1,1,100);
    [X Y]= meshgrid(x,y);
    subplot(2,2,1)
    Z= sin(n*pi*X).*sin(m*pi*Y);
    surf(X,Y,Z)
    colorbar
    shading interp
    
    subplot(2,2,2)
    Z= sin(n*pi*X).*cos((m+1/2)*pi*Y);
    surf(X,Y,Z)
    colorbar
    shading interp
    
    subplot(2,2,3)
    Z= cos((n+.5)*pi*X).*sin(m*pi*Y);
    surf(X,Y,Z)
    colorbar
    shading interp
    
    subplot(2,2,4)
    Z= cos((n+.5)*pi*X).*cos((m+.5)*pi*Y);
    surf(X,Y,Z)
    colorbar
    shading interp
    
return
    