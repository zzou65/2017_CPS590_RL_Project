function [] = run_disp_p(maxp,nvert)
%RUN_DISP_P run disp_p for various methods
%
%   RUN_DISP_P(MAXP,NVERT) runs DISP_P with default parameters for up to
%   MAXP local plane wave basis functions for elements with NVERT edges.
%
%   See also disp_p.

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
  omega_0 = [0.5 1 2 4 8];
  omega_1 = [0.125 0.25 0.5];
  
  % Run disp_p for ultra-weak variational method with omega scale 0
  flux_params = {};
  mtd_name = 'Ultra Weak Method';
  fig = disp_p(ndir,omega_0,0,flux_params,nvert,mtd_name);
  print(fig,'-depsc','disp_ultra_weak_0');
  saveas(fig,'disp_ultra_weak_0','fig');
  close(fig);
  
  % Run disp_p for ultra-weak variational method with omega scale 1
  flux_params = {};
  mtd_name = 'Ultra Weak Method';
  fig = disp_p(ndir,omega_1,1,flux_params,nvert,mtd_name);
  print(fig,'-depsc','disp_ultra_weak_1');
  saveas(fig,'disp_ultra_weak_1','fig');
  close(fig);
  
  % Run disp_p for PWDG with omega scale 0
  flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
  mtd_name = 'PWDG (\alpha=2/h\omega,\beta=0.01h\omega)';
  fig = disp_p(ndir,omega_0,0,flux_params,nvert,mtd_name);
  print(fig,'-depsc','disp_pwdg_0');
  saveas(fig,'disp_pwdg_0','fig');
  close(fig);
  
  % Run disp_p for PWDG with omega scale 1
  flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
  mtd_name = 'PWDG (\alpha=2/h\omega,\beta=0.01h\omega)';
  fig = disp_p(ndir,omega_1,1,flux_params,nvert,mtd_name);
  print(fig,'-depsc','disp_pwdg_1');
  saveas(fig,'disp_pwdg_1','fig');
  close(fig);
  
  % Reset path
  path(path0);
  
return