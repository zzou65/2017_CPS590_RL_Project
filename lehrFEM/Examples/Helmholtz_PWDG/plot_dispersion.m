function fig = plot_dispersion(ndir,omega,nVert,stretch,flux_params)
%PLOT_DISPERSION plot dispersion of PWDG
%
%   PLOT_DISPERSION(NDIR,OMEGA) plots the dispersion for various angles for
%   a discontinuous plane wave discretization with NDIR shape functions of 
%   the Helmholtz equation with wave number OMEGA.
%
%   PLOT_DISPERSION(NDIR,OMEGA,NVERT) uses a mesh of polygons with NVERT
%   edges for NVERT equal to three (default) or four.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Set default arguments
if(nargin<3 || isempty(nVert))
  nVert = 3;
end
if(nargin<4 || isempty(stretch))
  stretch = 1.1;
end
if(nargin<5 || isequal(flux_params,[]))
  flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
end

% Define discretization parameters
% params_in = {@(h,omega)1/(h*omega),@(h,omega)0.1*h*omega};
phi = 2*pi*(0:1/ndir:1-1/ndir)';

% Define mesh of angles
ntheta = 100;
theta = 2*pi*(0:1/ntheta:1-1/ntheta);

% Calculate disperison
k = dispersion(ndir,omega,theta,flux_params,nVert);

% Prepare for plotting
k(ntheta+1) = k(1);
theta = [theta,2*pi];
if(phi(1)==0)
  phi = [phi;2*pi];
end

% Plot wave number omega versus theta
fig = figure;
plot3(theta,real(k)/omega,imag(k)/omega,'b-');
hold on;
plot3(theta,ones(size(theta)),zeros(size(theta)),'r:');
plot3(phi,ones(size(phi)),zeros(size(phi)),'ko');
hold off;

% Set axis limits
set(gca,'XLim',[0,2*pi],'XTick',(0:4)/2*pi,'XTickLabel',{'0','1/2','1','3/2','2'});
rmax = max(abs(real(k)/omega-1));
imax = max(abs(imag(k)/omega));
set(gca,'YLim',1+stretch*rmax*[-1 1],'ZLim',stretch*imax*[-1,1]);

% Annotate plot
legend('num. sol.','exact','basis fn.');
xlabel('\bf \theta/\pi');
ylabel('\bf Re |k|/\omega');
zlabel('\bf Im |k|/\omega');
title(sprintf('\\bf Wave number k = |k|(cos(\\theta),sin(\\theta)) of PWDG solution of \\Delta u + \\omega^2 u=0, h\\omega=%g',omega));


% % Plot k versus theta
% figure;
% plot(theta,k/omega,'-',theta,ones(size(theta)),':',phi,ones(size(phi)),'ko');
% legend('numerical sol.','exact','basis fn.');
% set(gca,'XLim',[0,2*pi],'XTick',(0:4)/2*pi,'XTickLabel',{'0','1/2','1','3/2','2'});
% xlabel('\bf \theta/\pi');
% ylabel('\bf |k|/\omega');
% title(sprintf('\\bf Wave number k = |k|(cos(\\theta),sin(\\theta)) of PWDG solution of \\Delta u + \\omega^2 u=0, h\\omega=%g',omega));
% % title('\bf Dispersion Analysis of PWDG for \Delta u + \omega^2 u = 0');

% % Plot k vectors for which stiffness matrix is singular
% figure;
% plot(k.*cos(theta),k.*sin(theta),'-',...
%   omega*cos(theta),omega*sin(theta),':',...
%   omega*cos(phi),omega*sin(phi),'ko');
% legend('numerical sol.','exact','basis fn.');
% title('\bf Wave number k of PWDG solution of \Delta u + \omega^2 u = 0');
% 
% % Plot error in k
% figure
% semilogy(theta,abs(k-omega)/omega);
% set(gca,'XLim',[0,2*pi],'XTick',(0:4)/2*pi,'XTickLabel',{'0','1/2','1','3/2','2'});
% xlabel('\bf \theta/\pi');
% ylabel('\bf ||k|-\omega|/\omega');
% title('\bf Wave number k = |k|(cos(\theta),sin(\theta)) of PWDG solution of \Delta u + \omega^2 u = 0');
