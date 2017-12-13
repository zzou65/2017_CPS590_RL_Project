function [] = mg_eig_h()
% spectral radius of multigrid error propagation operator for various h
%
%   plots the spectral radius of the multigrid iteration for a sample
%   problem against the number of refinements of the finest mesh.  The
%   curve is bounded away from one, so the multigrid convergence rate is
%   independent of the mesh width.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data

  smoother = @gs_smooth;                         % multigrid smoother
  cyc = 1;                                       % V or W multigrid cycles
  m = [1 1];                                     % number of smoothing steps
  
  numits = 50;        % number of iterations in power method
  
  ref0 = 2;           % number of initial mesh refinements
  refs = 3:8;         % total number of mesh refinements
  
  % prealocate memory
  
  numrefs = length(refs);
  rho = nan(1,numrefs);
  
  % iterate over number of refinements
  
  for j=1:numrefs
    
    % generate multigrid data structure
    
    mg_data = mg_mesh('mesh',CMesh,'ref',[ref0,refs(j)]);
    mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
    mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother);
    mg_data = mg_error(mg_data,'iter',false);
    
    A = mg_data{end}.A;
    b = mg_data{end}.b;
    
    % calculate convergence rate for multigrid using the power method;
    % the error propagation operator is E = I - BA
    
    Eerr = ones(size(b));
    u0 = zeros(size(b));
    for i=1:numits
      err = Eerr;
      c = mg(u0,A*err,mg_data,0,1); % c = BAerr
      Eerr = err - c;
    end
    rho(j) = (err'*Eerr)/(err'*err);
    
  end
  
  % plot spectral radius vs. number of refinements
  
  figure;
  plot(refs,rho,'-d');
  grid on;
  set(gca,'XTick',refs);
  xlabel('\bf number of mesh refinements');
  ylabel('\bf maximal eigenvalue of error propagation matrix');
  title('\bf Dependence of Multigrid Convergence Rate on Mesh Size');
  
return