function [] = run_dispinv_p(maxp,nvert)
%RUN_DISPINV_P run dispinv_p
%
%   RUN_DISP_P(MAXP,NVERT) runs DISPINV_P with default parameters for up to
%   MAXP local plane wave basis functions for elements with NVERT edges.
%
%   See also dispinv_p.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Add parent directory to path
  path0 = path;
  cd ..;
  addpath(pwd);
  cd Experiments;

  % Define default arguments
  if(nargin<1 || isempty(maxp))
    maxp = 7;
  end
  if(nargin<2 || isempty(nvert))
    nvert = 3;
  end
  
  % Define parameters
  ndir = 3:maxp;
  omega = [0.5 1 2 4 8];
  
  % Run disp_p for ultra-weak variational method
  flux_params = {};
  mtd_name = 'Ultra Weak Method';
  fig = dispinv_p(ndir,omega,0,flux_params,nvert,mtd_name);
  print(fig,'-depsc','dispinv_ultra_weak');
  saveas(fig,'dispinv_ultra_weak','fig');
  close(fig);
  
  % Run disp_p for PWDG with omega scale 0
  flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
  mtd_name = 'PWDG (\alpha=2/h\omega,\beta=0.01h\omega)';
  fig = dispinv_p(ndir,omega,0,flux_params,nvert,mtd_name);
  print(fig,'-depsc','dispinv_pwdg');
  saveas(fig,'dispinv_pwdg','fig');
  close(fig);
  
  % Reset path
  path(path0);
  
return