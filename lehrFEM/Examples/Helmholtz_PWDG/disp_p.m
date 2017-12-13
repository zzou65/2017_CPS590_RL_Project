function fig = disp_p(ndir,omega0,scale_omega,flux_params,nVert,mtd_name)
%DISP_P plot dispersion vs. p
%   
%   DISP_P(NDIR,OMEGA0,SCALE_OMEGA,FLUX_PARAMS,NVERT,MTD_NAME) plots the
%   maximal relative error in the wave number of a discrete plane wave
%   against the number of local plane wave basis functions.
%
%   NDIR is the number of local plane waves.  It should be a sorted vector
%   of natural numbers.
%
%   OMEGA0 is the initial value of OMEGA.  It can be a vector or a scalar.
%   The actual value of OMEGA is (NDIR^SCALE_OMEGA)*OMEGA0.  The mesh width
%   is always one, so OMEGA is equal to H*OMEGA.
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
if(nargin<2 || isempty(omega0))
  omega0 = 1;
end
if(nargin<3 || isempty(scale_omega))
  scale_omega = 0;
end
if(nargin<4 || isempty(flux_params))
  flux_params = {};
end
if(nargin<5 || isempty(nVert))
  nVert = 3;
end
if(nargin<6 || isempty(mtd_name))
  mtd_name = 'PWDG';
end

% initialize
nndir = numel(ndir);
nomega = numel(omega0);
disp = nan(nomega,nndir);

% loop over values of ndir
for i=1:nndir
  
  % determine omega
  omega = ndir(i)^scale_omega*omega0;
  
  % loop over omega
  for j=1:nomega
    
    % determine discrete wave number
    k = dispersion(ndir(i),omega(j),[],flux_params,nVert);
    
    % calculate error
    disp(j,i) = abs(k-omega(j))/omega(j);
    
  end

end

% plot results
fig = figure;
lgd = cell(1,nomega);
if(scale_omega)
  omega_str = '\omega_0';
  scale_omega_str = rats(scale_omega);
  spaces = strfind(scale_omega_str,' ');
  scale_omega_str(spaces) = [];
  omega_title = sprintf(', %s = p^{%s}%s_0','\omega',scale_omega_str,'\omega');
else
  omega_str = '\omega';
  omega_title = '';
end
for j=1:nomega
  plot(ndir,disp(j,:));
  lgd{j} = sprintf('%s = %g',omega_str,omega0(j));
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
ylabel('\bf max_\theta |k(\theta) - \omega|/\omega');
title(sprintf('%s Dispersion of %s%s%s','\bf',mtd_name,lgd,omega_title));

