function [] = eig_linadv()
% plot spectral radius of multigrid methods for linear advection
%
%   PROBLEM : VALUES SEEM TO OSCILLATE !!

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % define parameters

  refs = [2,6];
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  
  f = @(x,varargin) zeros(size(x,1),1);
  gd = @(x,varargin) x(:,1)<0.5;
  v = @(x,varargin) ones(size(x,1),1)*[0.1,1];
  
  c = [1e-2,5e-3,1e-3,5e-4,1e-4];
  
  m = 1;
  
  maxit = 34;
  numc = length(c);
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.pre.its = m(1);
  AMGOpt.post.its = m(end);
  
  % initialize data
  
  eig_gmg = nan(1,numc);
  eig_gmgdw = nan(1,numc);
  eig_amg = nan(1,numc);
  eig_amgdw = nan(1,numc);
  
  % loop over values of diffusion parameter c
  
  for j=1:numc
    
    % generate geometric multigrid data structure
  
    mg_data = mg_mesh('mesh',CMesh,'ref',refs);
    mg_data = mg_stima(mg_data,'f',f,'gd',gd,'stima','assem',...
      'stima_assem',@(mesh) assem_stima(mesh,c(j),v));
    mg_data = mg_smooth(mg_data,'m',m,'smoother',@gs_smooth);
    mg_data = mg_error(mg_data,'energy',false,'iter',false);
    A = mg_data{end}.A;
    b = mg_data{end}.b;
    
    % calculate convergence rate for multigrid using the power method;
    % the error propagation operator is E = I - BA
    
    Eerr = ones(size(b));
    u0 = zeros(size(b));
    for i=1:maxit
      err = Eerr;
      corr = mg(u0,A*err,mg_data,0,1); % c = BAerr
      Eerr = err - corr;
      eig_gmg(j) = (err'*Eerr)/(err'*err);
    end
    
    % construct downwind smoother for geomtric multigrid
    
    v0 = v([0 0]);
    mg_data = mg_smooth(mg_data,'per','sort',...
      'per_fn',@(mesh,dofs) mesh.Coordinates(dofs,:)*v0');
    
    % calculate convergence rate for geom. multigrid with downwinding
    
    Eerr = ones(size(b));
    u0 = zeros(size(b));
    for i=1:maxit
      err = Eerr;
      corr = mg(u0,A*err,mg_data,0,1);
      Eerr = err - corr;
      eig_gmgdw(j) = (err'*Eerr)/(err'*err);
    end
    
    
    % generate algebraic multgrid data structure
    
    AMGData = AMGSetup(A,AMGOpt);
    
    % calculate convergence rate for algebraic multigrid using the power method

    Eerr = ones(size(b));
    for i=1:maxit
      err = Eerr;
      corr = AMGVcycle(AMGData,A*err);
      Eerr = err - corr;
      eig_amg(j) = (err'*Eerr)/(err'*err);
    end
    
    % construct downwind smoother for algebraic multigrid
    
    AMGOptDW = AMGOpt;
    [dummy,per] = sort(mg_data{end}.mesh.Coordinates(mg_data{end}.dofs,:)*v0');
    AMGOptDW.pre.GSperm = per;
    
    AMGDataDW = AMGSetup(A,AMGOptDW);
    
    % calculate convergence rate for algebraic multigrid using the power method

    Eerr = ones(size(b));
    for i=1:maxit
      err = Eerr;
      corr = AMGVcycle(AMGDataDW,A*err);
      Eerr = err - corr;
      eig_amgdw(j) = (err'*Eerr)/(err'*err);
    end
    
  end % for loop over values of diffusion parameter c, index j
  
  % plot eigenvalues against diffusion parameter
  
  figure;
  semilogx(c,abs(eig_gmg),'-o',...
           c,abs(eig_gmgdw),'-^',...
           c,abs(eig_amg),'-+',...
           c,abs(eig_amgdw),'-x');
  xlabel('\bf diffusion coefficient');
  ylabel('\bf maximal eigenvalue of error propagation matrix');
  title('\bf Convergence of Multigrid Methods for Linear Advection');
  legend('geom. multigrid','geom. multigrid, dw smoother',...
    'alg. multigrid','alg. multigrid, dw smoother',...
    'Location','NorthEast');
  grid on;
  set(gca,'XDir','reverse','XLim',[min(c),max(c)]);
  
return