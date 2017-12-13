function fig = dispinv_p(ndir,k0,scale_k,flux_params,nvert,mtd_name)
%DISPINV_P plot dispersion vs. p
%   
%   DISPINV_P(NDIR,K0,SCALE_K,FLUX_PARAMS,NVERT,MTD_NAME) plots the
%   maximal relative error in the wave number of the Helmholtz equation
%   satisfied by a discrete plane wave against the number of local plane
%   wave basis functions.
%
%   NDIR is the number of local plane waves.  It should be a sorted vector
%   of natural numbers.
%
%   K0 is the initial value of K.  It can be a vector or a scalar.
%   The actual value of K is (NDIR^SCALE_K)*K0.  The mesh width
%   is always one, so K is equal to H*K.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   NVERT is the number of vertices in an element of the mesh, ie. 3 for
%   a triangular mesh (default) and 4 for a square mesh.
%
%   MTD_NAME is a string containing the name of the method used, for
%   example 'PWDG'.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% read arguments
if(nargin<1 || isempty(ndir))
  ndir = 3:8;
end
if(nargin<2 || isempty(k0))
  k0 = 1;
end
if(nargin<3 || isempty(scale_k))
  scale_k = 0;
end
if(nargin<4 || isempty(flux_params))
  flux_params = {};
end
if(nargin<5 || isempty(nvert))
  nvert = 3;
end
if(nargin<6 || isempty(mtd_name))
  mtd_name = 'PWDG';
end

% initialize
nndir = numel(ndir);
nomega = numel(k0);
disp = nan(nomega,nndir);
% ntheta = 100;
% theta = 2*pi*(0:1/ntheta:1-1/ntheta)';
% k_dir = [cos(theta),sin(theta)];

% loop over values of ndir
for i=1:nndir
  
  % determine wave number
  k_ = ndir(i)^scale_k*k0;
  
  % construct wave directions between basis functions
  ntheta = ndir(i);
  theta = pi*(1/ntheta:2/ntheta:2-1/ntheta)';
  k_dir = [cos(theta),sin(theta)];
  
  % loop over wave number
  for j=1:nomega
    
    % determine discrete wave number
    omegas = dispinv(ndir(i),k_(j)*k_dir,flux_params,nvert);
    
    % calculate error
    disp(j,i) = max(abs(omegas-k_(j))/k_(j));
    
  end

end

% plot results
fig = figure;
lgd = cell(1,nomega);
if(scale_k)
  k_str = 'k_0';
  scale_k_str = rats(scale_k);
  spaces = strfind(scale_k_str,' ');
  scale_k_str(spaces) = [];
  k_title = sprintf(', k = p^{%s}k_0',scale_k_str);
else
  k_str = 'k';
  k_title = '';
end
for j=1:nomega
  plot(ndir,disp(j,:));
  lgd{j} = sprintf('%s = %g',k_str,k0(j));
  hold all;
end
hold off;
set(gca,'YScale','log','XTick',ndir);
if(nomega>1)
  legend(lgd{:},'Location','SouthWest');
  lgd = '';
else
  lgd = sprintf(', %s',lgd{1});
end
grid on;
xlabel('\bf p : number of plane wave shape functions');
ylabel('\bf max_\theta |k - \omega|/k');
title(sprintf('\\bf Dispersion of %s%s%s',mtd_name,lgd,k_title));
