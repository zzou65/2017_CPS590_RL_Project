function fig = plot_dispinv(ndir,k_,ntheta,nvert,nomega,stretch,flux_params)
%PLOT_DISPINV plot dispersion of PWDG
%   
%   PLOT_DISPINV(NDIR,K) plots the dispersion for various angles for
%   a discontinuous plane wave discretization of the Helholtz equation with
%   NDIR local shape functions.  The basis functions have wave numbers 
%   K = |k| and the Bloch waves have wave vector k.
%
%   The three-dimensional line plot has the angle between the propagation
%   direction of the wave and the x-axis on the x-axis.  The value of the
%   real part of the numerically computed Helmholtz wave number is on the
%   y-axis and the imaginary part is on the z-axis.  Both of these are
%   scaled by K.
%
%   FIG = PLOT_DISPINV(NDIR,K,NTHETA,NVERT,NOMEGA,STRETCH,FLUX_PARAMS) is 
%   the syntax including all the optional arguments.
%
%   NDIR is the number of local plane wave basis functions.
%
%   K the scalar wave number of the Bloch wave under consideration and also
%   of the discontinuous plane wave basis functions.
%
%   NTHETA is the number of propagation directions for which to determine
%   the dispersion relation.  The default value is 100.
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   NOMEGA is the number of wave numbers of the Helmholtz equation
%   satisfied by the Bloch wave to determine.  For small K, only one of
%   these is near K, but for large K, multiple values may exist.  The
%   default value is 1.  To plot possible values, an arbitrarily large
%   value is allowed, for example Inf.
%
%   The parameter STRETCH determines how tightly the axis limits fit the
%   plot (for the first Helmholtz wave number).
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are
%       FLUX_PARAMS = {'a',@(omega,h,varargin) 2/(h*omega),...
%         'b',@(omega,h,varargin) 0.01*h*omega};
%   which seems to produce relatively little dispersion.
%
%   The output argument FIG is an index pointer to the figure containing
%   the plot created by this function.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Read arguments
if(nargin<3 || isempty(ntheta))
  ntheta = 100;
end
if(nargin<4 || isempty(nvert))
  nvert = 3;
end
if(nargin<5 || isempty(nomega))
  nomega = 1;
end
if(nargin<6 || isempty(stretch))
  stretch = 1.1;
end
if(nargin<7 || isequal(flux_params,[]))
  flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
end

% Define discretization parameters
phi = 2*pi*(0:1/ndir:1-1/ndir)';

% Define mesh of angles
theta = 2*pi*(0:1/ntheta:1-1/ntheta)';
k = k_*[cos(theta),sin(theta)];

% Calculate disperison
omega = dispinv(ndir,k,flux_params,nvert,nomega);
nomega = size(omega,2);

% Prepare for plotting
omega(ntheta+1,:) = omega(1,:);
theta = [theta;2*pi];
if(phi(1)==0)
  phi = [phi;2*pi];
end

% Plot wave number omega versus theta
fig = figure;
for j=1:nomega
  h1 = plot3(theta,real(omega(:,j))/k_,imag(omega(:,j))/k_,...
    'Color','b','LineStyle','-','LineWidth',2.^(1-j));
  hold on;
end
h2 = plot3(theta,ones(size(theta)),zeros(size(theta)),'r:');
h3 = plot3(phi,ones(size(phi)),zeros(size(phi)),'ko');
hold off;

% Set axis limits
set(gca,'XLim',[0,2*pi],'XTick',(0:4)/2*pi,'XTickLabel',{'0','1/2','1','3/2','2'});
rmax = max(abs(real(omega(:,1))/k_-1));
imax = max(abs(imag(omega(:,1))/k_));
set(gca,'YLim',1+stretch*rmax*[-1 1],'ZLim',stretch*imax*[-1,1]);

% Annotate plot
legend([h1,h2,h3],'num. sol.','exact','basis fn.');
xlabel('\bf \theta/\pi');
ylabel('\bf Re \omega/|k|');
zlabel('\bf Im \omega/|k|');
title(sprintf('\\bf Wave number \\omega of Helmoltz eq. satisfied by k = |k|(cos(\\theta),sin(\\theta)) - periodic PWDG, |k|=%g',k_));

return