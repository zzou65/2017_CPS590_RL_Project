function [] = run_pollut_omega(maxrefs,max_omega_p,p,nvert)
%RUN_POLLUT_OMEGA run pollut_omega for various methods
%   
%   RUN_POLLUT_OMEGA(MAXREFS,MAX_OMEGA_P,P,NVERT) runs POLLUT_OMEGA with
%   1 to MAXREFS refinements of the initial mesh and omega equal to the
%   powers of two between 8 and 2^MAX_OMEGA_P.  P plane wave shape
%   functions are used (P can be a vector).  NVERT is the number of edges
%   of the elements in the mesh.
%
%   See also pollut_omega.

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
  if(nargin<1 || isempty(maxrefs))
    maxrefs = 4;
  end
  if(nargin<2 || isempty(max_omega_p))
    max_omega_p = 4;
  end
  if(nargin<3 || isempty(p))
    p = 4;
  end
  if(nargin<4 || isempty(nvert))
    nvert = 3;
  end
  
  % Define parameters
  refs = [1,maxrefs];
  omega = 2.^(3:max_omega_p);
  
  % Loop over p
  for j=1:numel(p)
    
    p_str = num2str(p(j));
  
    % Run pollut_omega for ultra-weak variational method
    flux_params = {};
    mtd_name = 'Ultra Weak Method';
    [fig1,fig2,fig3] = pollut_omega(refs,p(j),nvert,omega,[],flux_params,mtd_name);

    % Save figures
    print(fig1,'-depsc',['pollution_uwvf_maxerr_',p_str]);
    saveas(fig1,['pollution_uwvf_maxerr_',p_str],'fig');
    close(fig1);
    print(fig2,'-depsc',['pollution_uwvf_startconv_',p_str]);
    saveas(fig2,['pollution_uwvf_startconv_',p_str],'fig');
    close(fig2);
    print(fig3,'-depsc',['pollution_uwvf_plot_',p_str]);
    saveas(fig3,['pollution_uwvf_plot_',p_str],'fig');
    close(fig3);

    % Run pollut_omega for PWDG method
    flux_params = {'a',@(omega,h,varargin)1/(h*omega),...
      'b',@(omega,h,varargin)0.1*h*omega,'d',1};
    mtd_name = 'PWDG (\alpha=1/h\omega,\beta=0.1h\omega,\delta=1)';
    [fig1,fig2,fig3] = pollut_omega(refs,p(j),nvert,omega,[],flux_params,mtd_name);

    % Save figures
    print(fig1,'-depsc',['pollution_pwdg_maxerr_',p_str]);
    saveas(fig1,['pollution_pwdg_maxerr_',p_str],'fig');
    close(fig1);
    print(fig2,'-depsc',['pollution_pwdg_startconv_',p_str]);
    saveas(fig2,['pollution_pwdg_startconv_',p_str],'fig');
    close(fig2);
    print(fig3,'-depsc',['pollution_pwdg_plot_',p_str]);
    saveas(fig3,['pollution_pwdg_plot_',p_str],'fig');
    close(fig3);

  end
  
  % Reset path
  path(path0);

return