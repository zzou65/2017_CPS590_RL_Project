function [] = complexity_cascadic()
% complexity of cascadic multigrid
%
%   plots the time required by cascadic multigrid to solve a sample problem
%   vs. the number of degrees of freedom on the fines mesh for various
%   amounts of refinement and the least-squares line through these points.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % initialize constants
 
  ref = [1 8];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  du_ex = @(x,varargin) 2*x;                     % derivative of exact solution
  smoother = @sgs_smooth;                        % multigrid smoother
  tol0 = 5e-1;                                   % Stopping criterion
  cyc = 0;                                       % V or W multigrid cycles
  m = [1 0];                                     % number of smoothing steps
  M = 0.25;                                      % constant to determine maximal number of iterations per level
  Q = 2;                                         % constant to determine maximal number of iterations per level
  R = 5e-1;                                      % factor for tolerance on coarse grids

  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother);
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false,'ctrl','');
%   mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false); % use error estimator
  LVL = length(mg_data);
  LVL1 = LVL-4;
  
  % determine number of degrees of freedom on all levels
  
  dofs = zeros(1,LVL);
  for lvl=1:LVL
    dofs(lvl) = mg_data{lvl}.n.free;
  end
  
  % time cascadic multigrid for various levels of refinement
 
  numsgs = zeros(1,LVL-LVL1+1);
  disc_errs = zeros(1,LVL-LVL1+1);
  trunc_errs = zeros(1,LVL-LVL1+1);
  for lvl=LVL1:LVL
    
    % determine discretization error
    u = mg_data{lvl}.u_bd;
    u(mg_data{lvl}.dofs) = mg_data{lvl}.A\mg_data{lvl}.b;
    disc_err = H1SErr_LFE(mg_data{lvl}.mesh,u,P7O6(),du_ex);
    
    % deterimine maximal number of iterations on each level and tolerance
    maxit = ceil(M*lvl^2*Q.^(2*(lvl-(2:lvl))));
%     maxit = ceil(M*Q.^(2*(lvl-(2:lvl))));
%     tol = tol0*disc_err*[R*ones(1,lvl-2),1];
    tol = tol0*disc_err*[R*2.^(lvl-3:-1:0),1];
    
    m_tot = m(1)+m(end);
    maxit = ceil(maxit/m_tot);
    
    % calculate solution with cascadic multigrid

    [x,conv] = nmg_solve({mg_data{1:lvl}},tol,maxit);
    
    % count number of point sgs steps
    
    numsgs(lvl-LVL1+1) = 2*m_tot*conv.iter*dofs(2:lvl)';
    
    % store convergence information
    
    disc_errs(lvl-LVL1+1) = disc_err;
    trunc_errs(lvl-LVL1+1) = conv.error{end}.energy_exact(end);

    % determine truncation errors at the end of each iteration

    errs = zeros(1,lvl-1);
    for i=1:lvl-1
      errs(i) = conv.error{i}.energy_exact(end);
    end
    
    % display some convergence information

    fprintf(['Dofs       :',sprintf(' %10.0d ',dofs(2:lvl)),'\n']);
    fprintf(['Iterations :',sprintf(' %10.0d ',conv.iter*m_tot),'\n']);
    fprintf(['Error      :',sprintf(' %1.3e ',errs),'\n\n']);
    
  end
  
  % plot time against number of degrees of freedom
  
  figure;
  plot(dofs(LVL1:LVL),numsgs,'bo');
  hold on;
  set(gca,'XScale','log','YScale','log');
  grid('on');

  title('{\bf Complexity of Cascadic Multigrid}');
  xlabel('{\bf degrees of freedom}');
  ylabel('{\bf number of point Gauss-Seidel steps}')

  p = polyfit(log(dofs(LVL1:LVL)),log(numsgs),1);
  x = get(gca,'XLim');
  y = get(gca,'YLim');
  plot(x,exp(polyval(p,log(x))));
  set(gca,'YLim',y);
  add_Slope(gca,'NorthWest',p(1));
  
  % plot discretization and truncation errors
  
  figure;
  plot(dofs(LVL1:LVL),disc_errs,'rx-',dofs(LVL1:LVL),trunc_errs,'bo-');
  set(gca,'XScale','log','YScale','log');
  grid on;
  
  title('\bf Errors of Cascadic Multigrid Iteration');
  xlabel('\bf degrees of freedom');
  ylabel('\bf error in energy norm');
  legend('discretization error','truncation error','Location','NorthEast');
  
%   % plot truncation error over discretization error
%   
%   figure;
%   plot(dofs(LVL1:LVL),trunc_errs./disc_errs,'bs-');
%   set(gca,'XScale','log');
%   grid on;
%   
%   title('\bf Errors of Cascadic Multigrid Iteration');
%   xlabel('\bf degrees of freedom');
%   ylabel('\bf trunc. error over disc. error in energy norm');
  
return  