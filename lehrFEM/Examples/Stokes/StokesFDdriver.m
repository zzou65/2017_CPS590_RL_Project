function [ ] = StokesFDdriver()
% Driver function for naive finite differences for Stokes problem

    N = 50; % Griod resolution
    [u1,u2,p] = StokesFD(N,@f); % finite difference computation

    [X,Y] = meshgrid(0.05:0.05:0.95,0.05:0.05:0.95);
    fig = figure(); axis([-0.1 1.1 -0.1 1.1]); hold on;
    quiver(X,Y,zeros(size(X)),cos(pi*X),4); 
    plot([0 1 1 0 0],[0 0 1 1 0],'k-');
    xlabel('{\bf x_1}','fontsize',14);
    ylabel('{\bf x_2}','fontsize',14);
    print(fig, '-depsc2', sprintf('../../Slides/NPDEPics/stokesforce.eps', N) );

    unk = size(u1,1); h = 1/N; x = h:h:1-h; [X,Y] = meshgrid(x,x);
    fig = figure(); axis([-0.1 1.1 -0.1 1.1]); hold on;
    quiver(X,Y,u1,u2,4); plot([0 1 1 0 0],[0 0 1 1 0],'k-');
    xlabel('{\bf x_1}','fontsize',14);
    ylabel('{\bf x_2}','fontsize',14);
    print(fig, '-depsc2', sprintf('../../Slides/NPDEPics/stokesvres%d.eps', N) );

    [X,Y] = meshgrid(0:h:1,0:h:1);
    p = [zeros(1,unk+2);zeros(unk,1) , p , zeros(unk,1); zeros(1,unk+2)];
    fig2 = figure(); 
    mesh(X,Y,p);
    print(fig2, '-depsc2', sprintf('../../Slides/NPDEPics/stokespres%d.eps', N) );
end

function [force] = f(x1,x2)
    force = zeros(2,1); 
    force(2) = cos(pi*x1);
end
