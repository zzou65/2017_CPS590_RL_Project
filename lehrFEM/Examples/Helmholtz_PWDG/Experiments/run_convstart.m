function [] = run_convstart(maxrefs,lambda_c,ndir,d)
%RUN_CONVSTART run convstart for various methods
%   
%   RUN_CONVSTART(MAXREFS,LAMBDA_C,NDIR,D) runs the CONVSTART script for
%   three versions of the plane wave discontinuous Galerkin method and
%   saves the resulting figures in .EPS and .FIG formats.
%
%   MAXREFS is the maximal number of mesh refinements performed.  The
%   default value is 4.
%
%   LAMBDA_C is the frequency of the exact solution, that is, the exact
%   solution of the sample problem is a sine wave with angular frequency
%   LAMBDA_C times omega.  The default value is 1, which leads to an exact
%   solution in the kernel of the Helmholtz operator.
%
%   NDIR is the number of local plane wave basis functions; the default is
%   5.
%
%   D is the direction of propagation of the exact solution; the default
%   value is [2 1].
%
%   See also convstart.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Define default arguments
  if(nargin<1 || isempty(maxrefs))
    maxrefs = 4;
  end
  if(nargin<2 || isempty(lambda_c))
    lambda_c = 1;
  end
  if(nargin<3 || isempty(ndir))
    ndir = 4;
%     ndir = 5;
  end
  if(nargin<4 || isempty(d))
    d = [1 1];
%     d = [2 1];
  end
  
  % Define parameters
  refs = [1,maxrefs];
  omega_p = 3:7;
  lambda_ = @(omega) lambda_c*omega;
  
  % Run convstart for ultra-weak variational method
  flux_params = {};
  mtd_name = ['Ultra Weak Method (p=',num2str(ndir),')'];
  [fig1,fig2,homega,err_discr] = convstart(refs,omega_p,ndir,lambda_,d,...
    flux_params,mtd_name);
  
  % Save figures
  print(fig1,'-depsc','ultra_weak_1');
  saveas(fig1,'ultra_weak_1','fig');
  close(fig1);
  print(fig2,'-depsc','ultra_weak_2');
  saveas(fig2,'ultra_weak_2','fig');
  close(fig2);
  
  % Run convstart for first PWDG method
  flux_params = {'a',@(omega,h,varargin)1/(h*omega),'b',0,'d',1};
  mtd_name = ['PWDG (\alpha=1/h\omega,\beta=0,\delta=1,p=',num2str(ndir),')'];
  [fig1,fig2] = convstart(refs,omega_p,ndir,lambda_,d,...
    flux_params,mtd_name,homega,err_discr);
  
  % Save figures
  print(fig1,'-depsc','pwdg_a_1');
  saveas(fig1,'pwdg_a_1','fig');
  close(fig1);
  print(fig2,'-depsc','pwdg_a_2');
  saveas(fig2,'pwdg_a_2','fig');
  close(fig2);
  
  % Run convstart for second PWDG method
  flux_params = {'a',@(omega,h,varargin)1/(h*omega),...
    'b',@(omega,h,varargin)0.1*h*omega,'d',1};
  mtd_name = ['PWDG (\alpha=1/h\omega,\beta=0.1h\omega,\delta=1,p=',num2str(ndir),')'];
  [fig1,fig2] = convstart(refs,omega_p,ndir,lambda_,d,...
    flux_params,mtd_name,homega,err_discr);
  
  % Save figures
  print(fig1,'-depsc','pwdg_b_1');
  saveas(fig1,'pwdg_b_1','fig');
  close(fig1);
  print(fig2,'-depsc','pwdg_b_2');
  saveas(fig2,'pwdg_b_2','fig');
  close(fig2);
  
return